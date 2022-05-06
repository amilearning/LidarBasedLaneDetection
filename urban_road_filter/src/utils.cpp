#include "urban_road_filter/data_structures.hpp"


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


auto Detector::polyfit(std::vector<double> x, std::vector<double> y){
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
    return f.getCoefficients();
  }
  catch (std::exception &e)
  {
     ROS_WARN("PolyFit Error");    
    Eigen::Matrix<double, -1, 1> NA; 
    return NA;
   
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
        auto filterCondition = boost::make_shared<FilteringCondition<pcl::PointXYZI>>(
            [=](const pcl::PointXYZI& point){
                return point.x >= laneStartX[i] && laneStartX[i+1] &&
                point.y >= laneStartY - horizontalBinResolution/2 && point.y <= laneStartY + horizontalBinResolution/2;
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
               
              if(value_x.size() >= 3){
                auto f = polyfit(value_x, value_y);
                if( f.size() == 0){
                  ROS_WARN("No solution found from polyfit"); 
                  return; 
                }
                std::cout << "coeefienct"<< std::endl;
                std::cout << f << std::endl;
                std::cout << "~~~"<< std::endl;
              }else if(value_x.size() == 2){
                    ROS_WARN("liner");
              }else{
                  verticalBins[i][0][j] = verticalBins[numVerticalBins-1][0][j];
                  verticalBins[i][1][j] = verticalBins[numVerticalBins-1][1][j];
                  verticalBins[i][2][j] = verticalBins[numVerticalBins-1][2][j];
              }
              
                
                
                   
                
                // error =  mean(sqrt((polyval(P, value(:, 1)) - value(:, 2)).^2));
                // if error < 0.1
                //     xval = (roi(1) + roi(2))/2;
                //     yval = polyval(P, xval);
                    
                //     % Use error to regularize the value of predicted y
                //     yval = yval - error*abs(yval);
                //     startLanePoints(j) = yval;
                //     roi(3:4) = [yval - horizontalBinResolution, ...
                //         yval + horizontalBinResolution];
                //     % Update the lane point with the centre of the
                //     % predicted window
                //     tmpPc = select(ptCloud, findPointsInROI(ptCloud, roi));
                //     zmean = mean(tmpPc.Location(:, 3));
                //     verticalBins(i, :, j) = [xval, yval, zmean];
                //     if display
                //         m = helperDrawCuboidFromROI(roi, ptCloud);
                //         ax = plot(m);
                //         ax.FaceColor = 'green';
                //         ax.FaceAlpha = 1;
                //     end
                // else
                //     roi(3:4) = [startLanePoints(j) - horizontalBinResolution, ...
                //         startLanePoints(j) + horizontalBinResolution];
                //     if display
                //         m = helperDrawCuboidFromROI(roi, ptCloud);
                //         ax = plot(m);
                //         ax.FaceColor = 'green';
                //         ax.FaceAlpha = 1;
                //     end
                // end
                
        }
        
      }

    }

    // lane1 = lanes(:, :, 1);
    // lane2 = lanes(:, :, 2);
    // lane1(all(lane1 == 0, 2), :) = [];
    // lane2(all(lane2 == 0, 2), :) = [];
    // detectedLanePoints{1} = lane1;
    // detectedLanePoints{2} = lane2;
}



        
