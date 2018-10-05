#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <thread>
#include "Eigen-3.3/Eigen/Core"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }

double deg2rad(double x) { return x * pi() / 180; }

double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
    auto found_null = s.find("null");
    auto b1 = s.find_first_of("[");
    auto b2 = s.find_first_of("}");
    if (found_null != string::npos) {
        return "";
    } else if (b1 != string::npos && b2 != string::npos) {
        return s.substr(b1, b2 - b1 + 2);
    }
    return "";
}

double distance(double x1, double y1, double x2, double y2) {
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y) {

    double closestLen = 100000; //large number
    int closestWaypoint = 0;

    for (int i = 0; i < maps_x.size(); i++) {
        double map_x = maps_x[i];
        double map_y = maps_y[i];
        double dist = distance(x, y, map_x, map_y);
        if (dist < closestLen) {
            closestLen = dist;
            closestWaypoint = i;
        }

    }

    return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y) {

    int closestWaypoint = ClosestWaypoint(x, y, maps_x, maps_y);

    double map_x = maps_x[closestWaypoint];
    double map_y = maps_y[closestWaypoint];

    double heading = atan2((map_y - y), (map_x - x));

    double angle = fabs(theta - heading);
    angle = min(2 * pi() - angle, angle);

    if (angle > pi() / 4) {
        closestWaypoint++;
        if (closestWaypoint == maps_x.size()) {
            closestWaypoint = 0;
        }
    }

    return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y) {
    int next_wp = NextWaypoint(x, y, theta, maps_x, maps_y);

    int prev_wp;
    prev_wp = next_wp - 1;
    if (next_wp == 0) {
        prev_wp = maps_x.size() - 1;
    }

    double n_x = maps_x[next_wp] - maps_x[prev_wp];
    double n_y = maps_y[next_wp] - maps_y[prev_wp];
    double x_x = x - maps_x[prev_wp];
    double x_y = y - maps_y[prev_wp];

    // find the projection of x onto n
    double proj_norm = (x_x * n_x + x_y * n_y) / (n_x * n_x + n_y * n_y);
    double proj_x = proj_norm * n_x;
    double proj_y = proj_norm * n_y;

    double frenet_d = distance(x_x, x_y, proj_x, proj_y);

    //see if d value is positive or negative by comparing it to a center point

    double center_x = 1000 - maps_x[prev_wp];
    double center_y = 2000 - maps_y[prev_wp];
    double centerToPos = distance(center_x, center_y, x_x, x_y);
    double centerToRef = distance(center_x, center_y, proj_x, proj_y);

    if (centerToPos <= centerToRef) {
        frenet_d *= -1;
    }

    // calculate s value
    double frenet_s = 0;
    for (int i = 0; i < prev_wp; i++) {
        frenet_s += distance(maps_x[i], maps_y[i], maps_x[i + 1], maps_y[i + 1]);
    }

    frenet_s += distance(0, 0, proj_x, proj_y);

    return {frenet_s, frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double>
getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y) {
    int prev_wp = -1;

    while (s > maps_s[prev_wp + 1] && (prev_wp < (int) (maps_s.size() - 1))) {
        prev_wp++;
    }

    int wp2 = (prev_wp + 1) % maps_x.size();

    double heading = atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));
    // the x,y,s along the segment
    double seg_s = (s - maps_s[prev_wp]);

    double seg_x = maps_x[prev_wp] + seg_s * cos(heading);
    double seg_y = maps_y[prev_wp] + seg_s * sin(heading);

    double perp_heading = heading - pi() / 2;

    double x = seg_x + d * cos(perp_heading);
    double y = seg_y + d * sin(perp_heading);

    return {x, y};

}

int main() {
    uWS::Hub h;

    // Load up map values for waypoint's x,y,s and d normalized normal vectors
    vector<double> map_waypoints_x;
    vector<double> map_waypoints_y;
    vector<double> map_waypoints_s;
    vector<double> map_waypoints_dx;
    vector<double> map_waypoints_dy;

    // Waypoint map to read from
    string map_file_ = "../data/highway_map.csv";
    // The max s value before wrapping around the track back to 0
    double max_s = 6945.554;

    ifstream in_map_(map_file_.c_str(), ifstream::in);

    string line;
    while (getline(in_map_, line)) {
        istringstream iss(line);
        double x;
        double y;
        float s;
        float d_x;
        float d_y;
        iss >> x;
        iss >> y;
        iss >> s;
        iss >> d_x;
        iss >> d_y;
        map_waypoints_x.push_back(x);
        map_waypoints_y.push_back(y);
        map_waypoints_s.push_back(s);
        map_waypoints_dx.push_back(d_x);
        map_waypoints_dy.push_back(d_y);
    }

    // start in lane 1
    int lane = 1;

    // Have a reference velocity to target
    double refVel = 0.4; //mph

    // Keep a count of lane changes in the last 5 cycles
    int numLaneShifts = 1;

    h.onMessage(
            [&map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy, &lane, &numLaneShifts, &refVel](
                    uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                    uWS::OpCode opCode) {
                // "42" at the start of the message means there's a websocket message event.
                // The 4 signifies a websocket message
                // The 2 signifies a websocket event
                //auto sdata = string(data).substr(0, length);
                //cout << sdata << endl;
                if (length && length > 2 && data[0] == '4' && data[1] == '2') {

                    auto s = hasData(data);

                    if (s != "") {
                        auto j = json::parse(s);

                        string event = j[0].get<string>();

                        if (event == "telemetry") {
                            // j[1] is the data JSON object

                            // Main car's localization Data
                            double car_x = j[1]["x"];
                            double car_y = j[1]["y"];
                            double car_s = j[1]["s"];
                            double car_d = j[1]["d"];
                            double car_yaw = j[1]["yaw"];
                            double car_speed = j[1]["speed"];

                            // Previous path data given to the Planner
                            auto previous_path_x = j[1]["previous_path_x"];
                            auto previous_path_y = j[1]["previous_path_y"];
                            // Previous path's end s and d values
                            double end_path_s = j[1]["end_path_s"];
                            double end_path_d = j[1]["end_path_d"];

                            // Sensor Fusion Data, a list of all other cars on the same side of the road.
                            auto sensor_fusion = j[1]["sensor_fusion"];

                            int prev_size = previous_path_x.size();

                            json msgJson;

                            // Sensor Fusion
                            if (prev_size > 0) {
                                car_s = end_path_s;
                            }

                            // define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
                            bool tooClose = false;
                            bool leftLaneClear = true;
                            bool rightLaneClear = true;

                            // find ref_v to use
                            for (int i = 0; i < sensor_fusion.size(); i++) {
                                // car is in my lane
                                double vx = sensor_fusion[i][3];
                                double vy = sensor_fusion[i][4];
                                double checkSpeed = sqrt(vx * vx + vy * vy);
                                double checkCarS = sensor_fusion[i][5];

                                checkCarS += ((double) prev_size * 0.02 * checkSpeed);

                                double minFrontDist = 30.0; // safe distance from the front car
                                double minRearDist = 40.0; // safe distance from the rear car

                                // check car lane
                                float d = sensor_fusion[i][6];
                                if (d < (2 + 4 * lane + 2) && d > (2 + 4 * lane - 2)) {
                                    //check s values greater than min and s gap
                                    if ((checkCarS > car_s) && ((checkCarS - car_s) < minFrontDist)) {
                                        tooClose = true;
                                        if ((checkCarS > car_s) && ((checkCarS - car_s) < (minFrontDist * 0.6))) {
                                            refVel -= 0.70; // brake harder as it is close
                                            std::cout << "Hard Brake: Unsafe Distance" << std::endl;
                                        }
                                    }
                                }

                                if (tooClose) {
                                    // check for lane change possibility
                                    int rightLane = lane + 1;
                                    if (lane < 2 && rightLaneClear && d < (2 + 4 * rightLane + 2) &&
                                        d > (2 + 4 * rightLane - 2)) {
                                        if ((checkCarS > car_s) && ((checkCarS - car_s) < (minFrontDist * 1.3)) ||
                                            (checkCarS < car_s) && ((car_s - checkCarS) < minRearDist)) {
                                            rightLaneClear = false;
                                        }
                                    }

                                    int leftLane = lane - 1;
                                    if (lane > 0 && leftLaneClear && d < (2 + 4 * leftLane + 2) &&
                                        d > (2 + 4 * leftLane - 2)) {
                                        if ((checkCarS > car_s) && ((checkCarS - car_s) < (minFrontDist * 1.3)) ||
                                            (checkCarS < car_s) && ((car_s - checkCarS) < minRearDist)) {
                                            leftLaneClear = false;
                                        }
                                    }
                                }
                                numLaneShifts += 1;
                            }

                            if (tooClose) {
                                refVel -= 0.35;

                                if (lane > 0 && leftLaneClear && (numLaneShifts % 5 == 0)) {
                                    std::cout << "Steering Left" << std::endl;
                                    lane -= 1;
                                    numLaneShifts += 1;
                                } else if (lane < 2 && rightLaneClear && (numLaneShifts % 5 == 0)) {
                                    std::cout << "Steering Right" << std::endl;
                                    lane += 1;
                                    numLaneShifts += 1;
                                }

                            } else if (refVel <= 49.2) {
                                refVel += 0.5;
                            }



                            // List of widely spaced (x,y) waypoints, evenly spaced at 30m
                            vector<double> ptsx;
                            vector<double> ptsy;

                            // reference x,y,yaw states
                            double refX = car_x;
                            double refY = car_y;
                            double refYaw = deg2rad(car_yaw);

                            if (prev_size < 2) {
                                // make the path tangent to the car
                                double prevCarX = car_x - cos(car_yaw);
                                double prevCarY = car_y - sin(car_yaw);

                                ptsx.push_back(prevCarX);
                                ptsx.push_back(car_x);

                                ptsy.push_back(prevCarY);
                                ptsy.push_back(car_y);
                            } else {
                                // define reference states using previous path's end points
                                refX = previous_path_x[prev_size - 1];
                                refY = previous_path_y[prev_size - 1];

                                double refXPrev = previous_path_x[prev_size - 2];
                                double refYPrev = previous_path_y[prev_size - 2];
                                refYaw = atan2(refY - refYPrev, refX - refXPrev);

                                // make the path tangent to the previous path's end point
                                ptsx.push_back(refXPrev);
                                ptsx.push_back(refX);

                                ptsy.push_back(refYPrev);
                                ptsy.push_back(refY);
                            }


                            // Define a path that the car will visit every .02 seconds
                            // Add evenly spaced 30m points
                            vector<double> nextWayPoint0 = getXY(car_s + 30, (2 + 4 * lane), map_waypoints_s,
                                                            map_waypoints_x, map_waypoints_y);
                            vector<double> nextWayPoint1 = getXY(car_s + 60, (2 + 4 * lane), map_waypoints_s,
                                                            map_waypoints_x, map_waypoints_y);
                            vector<double> nextWayPoint2 = getXY(car_s + 90, (2 + 4 * lane), map_waypoints_s,
                                                            map_waypoints_x, map_waypoints_y);

                            ptsx.push_back(nextWayPoint0[0]);
                            ptsx.push_back(nextWayPoint1[0]);
                            ptsx.push_back(nextWayPoint2[0]);

                            ptsy.push_back(nextWayPoint0[1]);
                            ptsy.push_back(nextWayPoint1[1]);
                            ptsy.push_back(nextWayPoint2[1]);

                            for (int i = 0; i < ptsx.size(); i++) {
                                // reset car reference angle to 0deg
                                double shiftX = ptsx[i] - refX;
                                double shiftY = ptsy[i] - refY;

                                ptsx[i] = (shiftX * cos(0 - refYaw) - shiftY * sin(0 - refYaw));
                                ptsy[i] = (shiftX * sin(0 - refYaw) + shiftY * cos(0 - refYaw));
                            }

                            // Spline Creation
                            tk::spline s;

                            // set points for spline
                            s.set_points(ptsx, ptsy);

                            //define points to be used for planner
                            vector<double> nextXVals;
                            vector<double> nextYVals;

                            // resuse previous path points
                            for (int i = 0; i < previous_path_x.size(); i++) {
                                nextXVals.push_back(previous_path_x[i]);
                                nextYVals.push_back(previous_path_y[i]);
                            }

                            // select points on spline to travel at refVal
                            double targetX = 30.0;
                            double targetY = s(targetX);
                            double targetDist = sqrt((targetX) * (targetX) + (targetY) * (targetY));
                            double xAdditional = 0;

                            // Output 50 points
                            for (int i = 1; i <= 50 - previous_path_x.size(); i++) {
                                double N = (targetDist / (0.02 * refVel / 2.24));
                                double x_point = xAdditional + (targetX) / N;
                                double y_point = s(x_point);

                                xAdditional = x_point;

                                double xRef = x_point;
                                double yRef = y_point;

                                // reconvert ref. angle
                                x_point = (xRef * cos(refYaw) - yRef * sin(refYaw));
                                y_point = (xRef * sin(refYaw) + yRef * cos(refYaw));

                                x_point += refX;
                                y_point += refY;

                                nextXVals.push_back(x_point);
                                nextYVals.push_back(y_point);
                            }
                            msgJson["next_x"] = nextXVals;
                            msgJson["next_y"] = nextYVals;

                            auto msg = "42[\"control\"," + msgJson.dump() + "]";

                            ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

                        }
                    } else {
                        // Manual driving
                        std::string msg = "42[\"manual\",{}]";
                        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
                    }
                }
            });

    // We don't need this since we're not using HTTP but if it's removed the
    // program
    // doesn't compile :-(
    h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                       size_t, size_t) {
        const std::string s = "<h1>Hello world!</h1>";
        if (req.getUrl().valueLength == 1) {
            res->end(s.data(), s.length());
        } else {
            // i guess this should be done more gracefully?
            res->end(nullptr, 0);
        }
    });

    h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
        std::cout << "Connected!!!" << std::endl;
    });

    h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                           char *message, size_t length) {
        ws.close();
        std::cout << "Disconnected" << std::endl;
    });

    int port = 4567;
    if (h.listen(port)) {
        std::cout << "Listening to port " << port << std::endl;
    } else {
        std::cerr << "Failed to listen to port" << std::endl;
        return -1;
    }
    h.run();
}
