#include "lidar_lane_detector.h"

float       params::min_X,
            params::max_X,
            params::min_Y,
            params::max_Y,
            params::min_Z,
            params::max_Z;   

            
void marker_init(visualization_msgs::Marker& m)
{
    m.pose.position.x = 0;
    m.pose.position.y = 0;
    m.pose.position.z = 0;

    m.pose.orientation.x = 0.0;
    m.pose.orientation.y = 0.0;
    m.pose.orientation.z = 0.0;
    m.pose.orientation.w = 1.0;

    m.scale.x = 0.5;
    m.scale.y = 0.5;
    m.scale.z = 0.5;
}

inline std_msgs::ColorRGBA setcolor(float r, float g, float b, float a)
{
    std_msgs::ColorRGBA c;
    c.r = r;
    c.g = g;
    c.b = b;
    c.a = a;
    return c;
}

LaneDetector::LaneDetector(ros::NodeHandle* nh){
    /*subscribing to the given topic*/
    std::string point_cloud_topic_name; 
    nh->param<std::string>("pointcloud_topic", point_cloud_topic_name, "/os_cloud_node/points");
    sub = nh->subscribe(point_cloud_topic_name, 1, &LaneDetector::filtered,this);
    /*publishing filtered points*/
    pub_road = nh->advertise<pcl::PCLPointCloud2>("road", 1);
    pub_high = nh->advertise<pcl::PCLPointCloud2>("curb", 1);
    // pub_box = nh->advertise<pcl::PCLPointCloud2>("roi", 1); // ROI - region of interest
    pub_pobroad = nh->advertise<pcl::PCLPointCloud2>("road_probably", 1);
    pub_marker = nh->advertise<visualization_msgs::MarkerArray>("road_marker", 1);
    pub_point_debugger = nh->advertise<pcl::PCLPointCloud2>("roi", 1); // ROI - region of interest
    // LaneDetector::beam_init();

    ROS_INFO("Ready");

}

/*FUNCTIONS*/

void LaneDetector::filtered(const pcl::PointCloud<pcl::PointXYZI> &cloud){
    /*variables for the "for" loops*/
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    int i, j, k, l;
    
    pcl::PointXYZI pt;                                                                      //temporary variable for storing a point
    auto cloud_filtered_Box = boost::make_shared<pcl::PointCloud<pcl::PointXYZI>>(cloud);   //all points in the detection area
    pcl::PointCloud<pcl::PointXYZI> cloud_filtered_ROI;                                    //filtered points (ROI)
    pcl::PointCloud<pcl::PointXYZI> cloud_filtered_ProbablyLane;                         //filtered points (lanes)
    // pcl::PointCloud<pcl::PointXYZI> cloud_ground;                         //filtered points (ground)

       
    pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_ground(new pcl::PointCloud<pcl::PointXYZI>);


    auto filterCondition = boost::make_shared<FilteringCondition<pcl::PointXYZI>>(
        [=](const pcl::PointXYZI& point){
            return point.x >= params::min_X && point.x <= params::max_X &&
            point.y >= params::min_Y && point.y <= params::max_Y &&
            point.z >= params::min_Z && point.z <= params::max_Z &&
            point.x + point.y + point.z != 0;
        }
    );
    pcl::ConditionalRemoval<pcl::PointXYZI> condition_removal;
    condition_removal.setCondition(filterCondition);
    condition_removal.setInputCloud(cloud_filtered_Box);
    condition_removal.filter(*cloud_filtered_Box);

    /*number of points in the detection area*/
    size_t piece = cloud_filtered_Box->points.size();

    /*A minimum of 30 points are requested in the detection area to avoid errors.
    Also, there is not much point in evaluating less data than that.*/
    if (piece < 30){
        return;
    }
    
    pub_road.publish(*cloud_filtered_Box);
    // pub_point_debugger.publish(*cloud_filtered_Box);
    // ROS_INFO("cloud pub");
  

    //////////////////

     pcl::ModelCoefficients::Ptr coefficients (new pcl::ModelCoefficients ());
     pcl::PointIndices::Ptr inliers (new pcl::PointIndices ());


//   Create the segmentation object.
  pcl::SACSegmentation<pcl::PointXYZI> seg;
  seg.setOptimizeCoefficients (true);       //// Enable model coefficient refinement (optional).
  seg.setInputCloud (cloud_filtered_Box);                 // 
  seg.setModelType (pcl::SACMODEL_PLANE);    //// Configure the object to look for a plane.
  seg.setMethodType (pcl::SAC_MLESAC);       //  // Use MLESAC method.
  seg.setMaxIterations (1000);               
  seg.setDistanceThreshold (0.02);          // Set the maximum allowed distance to the model.
  //seg.setRadiusLimits(0, 0.1);     // For cylinder, Set minimum and maximum radii of the cylinder.
  seg.segment (*inliers, *coefficients);    // 


   

  // (eg. ax + by + cz + d = 0 ).
  std::cerr << "Model coefficients: " << coefficients->values[0] << " " 
                                      << coefficients->values[1] << " "
                                      << coefficients->values[2] << " " 
                                      << coefficients->values[3] << std::endl;

  pcl::copyPointCloud<pcl::PointXYZI>(*cloud_filtered_Box, *inliers, *cloud_ground);


    
    pub_point_debugger.publish(*cloud_ground);
    

std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
  std::cout << "loop iteration time = " << sec.count() << " seconds" << std::endl;

}