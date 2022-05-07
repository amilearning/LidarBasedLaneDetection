#include "urban_road_filter/data_structures.hpp"

void setColor(std_msgs::ColorRGBA* cl, double r, double g, double b, double a)
{
  cl->r = r;
  cl->g = g;
  cl->b = b;
  cl->a = a;
}


template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{
  
  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}


PolyFit<double> Detector::polyfit(std::vector<double> x, std::vector<double> y){
  try 
  {
    
    // Create options for fitting
    PolyFit<double>::Options options;
    options.polyDeg = 2; 
    options.solver = PolyFit<double>::EIGEN_JACOBI_SVD; // SVD solver
    options.maxMPts = x.size(); // x.size(); // #Pts to use for fitting
    options.maxTrial = 500; // Max. nr. of trials
    options.tolerance = params::Polyfit_tolerance; // Distance tolerance for inlier    
    // Create fitting object
    PolyFit<double>f(x,y,options);    
    // Solve using RANSAC 
    f.solveRLS(); 
     
    // Output result    
    polyfit_error = false;
    return f;
  }
  catch (std::exception &e)
  {
     ROS_WARN("PolyFit Error");  
     polyfit_error = true;  
    poly_error = 1e10;
  }
}

void Detector::computeHistogram(int numBins, pcl::PointCloud<pcl::PointXYZI>::Ptr cloudPtr, std::vector<double>& histval, std::vector<double>& yvals){

std::fill(histval.begin(), histval.end(), 0);
std::fill(yvals.begin(), yvals.end(), 0);

std::vector<double> binStartY = linspace(params::min_Y,params::max_Y,numBins);
    for (int i=0; i < numBins-1; i++){
        
        pcl::PassThrough<pcl::PointXYZI> pass;
        pcl::PointCloud<pcl::PointXYZI>::Ptr roiCloudPtr (new pcl::PointCloud<pcl::PointXYZI>);
        pass.setInputCloud (cloudPtr);           
        pass.setFilterFieldName ("y");  
        pass.setFilterLimits (binStartY[i], binStartY[i+1]);   
        pass.filter (*roiCloudPtr);                   
        if ( roiCloudPtr->points.size() > 0 ){
            double histval_tmp = 0.0;
            for (pcl::PointCloud<pcl::PointXYZI>::const_iterator it = roiCloudPtr->begin(); it != roiCloudPtr->end(); it++) {
                histval_tmp += it->intensity;            
            }
            histval[i] = histval_tmp;
            yvals[i] = (binStartY[i] + binStartY[i+1])/2.0;
        }

    }   


}





void Detector::lanewidthFilter(std::vector<double>& lane_filterd_yvals, std::vector<double>& detected_peak, std::vector<double> yvals,std::vector<double> peaks){
// yvals = startYs
// pkHistVal = peaks
  int min_i,min_j;
  int min_tmp = 1000; 
  double estimated_lane_width = -1.0; 

  std::vector<std::vector<int>> diff;
  for (size_t i = 0; i < yvals.size(); ++i ){ 
    if(yvals[i] >= 0){
      continue; // filter only left lanes
    }
    std::vector<int> row_;
    for (size_t j = 0; j < yvals.size(); ++j ){ 
      if(yvals[j] < 0){
        continue; // filter only right lanes
      }     
      double tmp_val = fabs(lanewidth - (yvals[i]- yvals[j]));
      row_.push_back(tmp_val);
      // diff[i][j] = fabs(lanewidth - (yvals[i]- yvals[j]));
      if(min_tmp > tmp_val){
        min_tmp = tmp_val;
        min_i = i;
        min_j = j;
        estimated_lane_width = tmp_val;
      }
    }
    diff.push_back(row_);
  }

  detected_peak.clear();
  lane_filterd_yvals.clear();

  if(estimated_lane_width < 0){      
    if(yvals.back() <= 0){ 
      // a lane is in left side --> extract the most left lane 
      lane_filterd_yvals.push_back(yvals.front()); 
      detected_peak.push_back(peaks.front());
    }else{
      // a lane is in right side --> extract the most right lane 
      lane_filterd_yvals.push_back(yvals.back()); 
      detected_peak.push_back(peaks.back());
    }
  }else{
    // left and right lanes are avaialble. 
    lane_filterd_yvals.push_back(yvals[min_i]);    
    lane_filterd_yvals.push_back(yvals[min_j]);
    detected_peak.push_back(peaks[min_i]);
    detected_peak.push_back(peaks[min_j]);
  }




}


void Detector::DetectLanes(pcl::PointCloud<pcl::PointXYZI>::Ptr cloudPtr,std::vector<double> startLanePoints){
    //  double horizontalBinResolution;
    //     double verticalBinResolution
    int numVerticalBins = ceil((params::max_X - params::min_X)/verticalBinResolution);
    int numLanes = startLanePoints.size();
    // 
    double verticalBins[numVerticalBins][3][numLanes];
    double lanes[numVerticalBins][3][numLanes];

      for (int i = 0; i < numVerticalBins; i++)
        for (int j = 0; j < 3; j++)
          for (int k=0 ; k < numLanes; k++){
            verticalBins[i][j][k] = 0.0;
            lanes[i][j][k] = 0.0;
          }
          

    std::vector<double> laneStartX = linspace(params::min_X,params::max_X,numVerticalBins);

    for(int i=0; i< numVerticalBins-1;i++){
      for (int j=0; j < numLanes; j++){
        double laneStartY = startLanePoints[j]; 
        // Define a vertical roi window 
        double roi[4] = {laneStartX[i], laneStartX[i+1],laneStartY - horizontalBinResolution/2,laneStartY + horizontalBinResolution/2};
    
        auto filterCondition = boost::make_shared<FilteringCondition<pcl::PointXYZI>>(
            [=](const pcl::PointXYZI& point){
                return point.x >= roi[0] && point.x <= roi[1] &&
                point.y >= roi[2] && point.y <= roi[3];
            }
        );
        pcl::PointCloud<pcl::PointXYZI>::Ptr roi_points_tmp(new pcl::PointCloud<pcl::PointXYZI>);
        pcl::ConditionalRemoval<pcl::PointXYZI> condition_removal;
        condition_removal.setCondition(filterCondition);
        condition_removal.setInputCloud(cloudPtr);
        condition_removal.filter(*roi_points_tmp);

        if(roi_points_tmp->points.size() > 0){
          // ROI selected  point cloud is not empty
              int maxIndex = 0;          
              int counter = 0;
              double max_intensity_tmp = 0.0;  
              pcl::PointXYZI max_point_location;        
              for (pcl::PointCloud<pcl::PointXYZI>::const_iterator it = roi_points_tmp->begin(); it != roi_points_tmp->end(); it++) {
                    if(max_intensity_tmp < it->intensity){
                      max_intensity_tmp = it->intensity;
                      maxIndex = counter;
                      max_point_location = *it;
                    }                
                    counter ++;           
                }
                verticalBins[i][0][j] = max_point_location.x;
                verticalBins[i][1][j] = max_point_location.y;
                verticalBins[i][2][j] = max_point_location.z;
                lanes[i][0][j] = max_point_location.x;
                lanes[i][1][j] = max_point_location.y;
                lanes[i][2][j] = max_point_location.z;    
                startLanePoints[j] = max_point_location.y; 
        }else{
         // For dash lanes update the sliding window by 2D polynomial              
              std::vector<double> value_x, value_y;              
              for(int vert_i = 0; vert_i < numVerticalBins-1;vert_i++){
                  if( fabs(lanes[vert_i+1][1][j]) > 1e-3){
                    std::vector<double> value_row_tmp; 
                    value_x.push_back(lanes[vert_i+1][0][j]);
                    value_y.push_back(lanes[vert_i+1][1][j]);                    
                  }         
                }
               
              if(value_x.size() >= 2){
                PolyFit<double> f = polyfit(value_x, value_y);
                if(polyfit_error){
                  ROS_WARN("No solution found from polyfit"); 
                  return; 
                }           
                auto tt = f.getCoefficients();
                poly_error = f.getMSE();         
                // std::cout << "coeefienct"<< std::endl;
                // std::cout << f << std::endl;
                // std::cout << "~~~"<< std::endl;
                if(poly_error < 0.1){

                    std::vector<double> xval;
                    xval.push_back((roi[0] + roi[1])/2);
                    std::vector<double> yval = f.evalPoly(xval);                    
                    // Use error to regularize the value of predicted y
                    yval[0] = yval[0] - poly_error*fabs(yval[0]);
                    startLanePoints[j] = yval[0];
                    roi[2] = yval[0] - horizontalBinResolution;
                    roi[3] = yval[0] + horizontalBinResolution;                    
                    // % Update the lane point with the centre of the
                    // % predicted window
                    auto filterCondition_update = boost::make_shared<FilteringCondition<pcl::PointXYZI>>(
                              [=](const pcl::PointXYZI& point){
                              return point.x >= roi[0] && point.x <= roi[1] &&
                              point.y >= roi[2] && point.y <= roi[3];
                              }
                        );                    
                    pcl::ConditionalRemoval<pcl::PointXYZI> condition_removal_update;
                    condition_removal_update.setCondition(filterCondition_update);
                    condition_removal_update.setInputCloud(cloudPtr);
                    condition_removal_update.filter(*roi_points_tmp);

                    double zmean;
                    int counter_update = 0;
                    for (pcl::PointCloud<pcl::PointXYZI>::const_iterator it = roi_points_tmp->begin(); it != roi_points_tmp->end(); it++) {
                          zmean+=it->intensity;          
                          counter_update ++;           
                      }
                    zmean = zmean / (counter_update+1e-6);
                    verticalBins[i][0][j] = xval[0];
                    verticalBins[i][1][j] = yval[0];
                    verticalBins[i][2][j] = zmean;
                }else{
                    roi[2] = startLanePoints[j] - horizontalBinResolution;
                    roi[3] = startLanePoints[j] + horizontalBinResolution;
                }
                  
              // }else if(value_x.size() == 2){
              //       ROS_WARN("liner");
              }else{
                  verticalBins[i][0][j] = verticalBins[numVerticalBins-1][0][j];
                  verticalBins[i][1][j] = verticalBins[numVerticalBins-1][1][j];
                  verticalBins[i][2][j] = verticalBins[numVerticalBins-1][2][j];
                  continue;
              }
                
        }
        
      }

    }
    
    
  std::vector<std::vector<double>> lanes_points_x, lanes_points_y, lanes_points_z;    
    
    for(int i=0 ;i <numLanes; i++){
      std::vector<double> lane_points_x, lane_points_y, lane_points_z; 
      for(int j=0; j< numVerticalBins; j++){
        if(fabs(lanes[j][0][i]) < 1e-4){
          continue;
        }
            lane_points_x.push_back(lanes[j][0][i]);
            lane_points_y.push_back(lanes[j][1][i]);
            lane_points_z.push_back(lanes[j][2][i]);                  
      }
      // check if it has enough number of points otherwise, do not add it as lanes 
      if(lane_points_x.size() > 2){
          lanes_points_x.push_back(lane_points_x);
          lanes_points_y.push_back(lane_points_y);
          lanes_points_z.push_back(lane_points_z);
      }
    }
  
   // Polynomial fitting of the detected lane points
    std::vector<std::vector<double>> final_lanes_points_x, final_lanes_points_y, final_lanes_points_z;    
    std::vector<double> xeval = linspace(params::min_X,params::max_X,80);
    std::vector<PolyFit<double>> fs;
    std::vector<double> poly_errors;
    double min_error;
    int min_idx;
    for(int i=0; i < lanes_points_x.size(); i++){
      PolyFit<double> tmpf = polyfit(lanes_points_x[i], lanes_points_y[i]);
      
                if(polyfit_error){
                  ROS_WARN("No solution found from polyfit"); 
                  return; 
                }          
      fs.push_back(tmpf);
      poly_errors.push_back(tmpf.getMSE());
      if(min_error > poly_errors.back()){
        min_error = poly_errors.back();
        min_idx = i;
      }
      std::vector<double> yeval = fs[i].evalPoly(xeval);
      final_lanes_points_x.push_back(xeval);
      final_lanes_points_y.push_back(yeval);
      std::vector<double> z_tmp(xeval.size(),0);
      final_lanes_points_z.push_back(z_tmp);
    }

    
    
    
   

  visualize_lanes(final_lanes_points_x,final_lanes_points_y,final_lanes_points_z);
  

  
}



void Detector::visualize_lanes(std::vector<std::vector<double>> lanes_points_x, std::vector<std::vector<double>> lanes_points_y,std::vector<std::vector<double>> lanes_points_z)
{

    visualization_msgs::MarkerArray lane_strips;
    std_msgs::ColorRGBA lane_color;
    setColor(&lane_color, 0.0, 1.0, 0.0, 1.0);
    

    visualization_msgs::Marker lane_strip;
    
    lane_strip.header.stamp = ros::Time::now();
    lane_strip.header.frame_id = params::fixedFrame;
    lane_strip.color = lane_color;
    lane_strip.action = visualization_msgs::Marker::ADD;          
    lane_strip.ns = "lanes";
    lane_strip.type = visualization_msgs::Marker::POINTS;
    lane_strip.scale.x = 1.0;
    lane_strip.scale.y = 1.0;
    lane_strip.scale.z = 1.0;
    // lane_strip.lifetime = ros::Duration(1.0);
    lane_strip.id = 200;
  // fill out lane line
  for (int i=0; i < lanes_points_x.size(); i++){
    for(int j=0; j < lanes_points_x[i].size(); j++){
      geometry_msgs::Point p;
    p.x = lanes_points_x[i][j];
    p.y = lanes_points_y[i][j];
    p.z = lanes_points_z[i][j];    
    
    lane_strip.points.push_back(p);
    }    
  }

  lane_strips.markers.push_back(lane_strip);

  pub_lanes_marker.publish(lane_strips);
}

