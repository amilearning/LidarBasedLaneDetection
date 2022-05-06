#include "lidar_lane_detector.h"

/*parameter setup*/
void paramsCallback(lidar_lane_detector::LidarFiltersConfig &config, uint32_t level){
    bool t = config.ground_filter_enable;
    double tt = config.ground_height;
    int ttt = config.curb_points;   
    params::min_X = config.min_x;
    params::max_X = config.max_x;
    params::min_Y = config.min_y;
    params::max_Y = config.max_y;
    params::min_Z = config.min_z;
    params::max_Z = config.max_z;
    ROS_INFO("Updated params %s", ros::this_node::getName().c_str());
}

/*MAIN*/
int main(int argc, char **argv)
{

    /*initializing ROS*/
    ros::init(argc, argv, "lidar_lane_detector");
    ROS_INFO("Initializing %s", ros::this_node::getName().c_str());

    /*code needed for GUI*/
    dynamic_reconfigure::Server<lidar_lane_detector::LidarFiltersConfig> server;
    dynamic_reconfigure::Server<lidar_lane_detector::LidarFiltersConfig>::CallbackType f;
    f = boost::bind(&paramsCallback, _1, _2);
    server.setCallback(f);

    /*NodeHandle*/
    ros::NodeHandle nh;
    LaneDetector lane_detector(&nh);

    ros::spin();
    return 0;
}