#include "ros/ros.h"
#include <math.h>
#include <vector>
#include <iostream>
#include <string>
#include <nav_msgs/Odometry.h>


/**
 * This tutorial demonstrates simple receipt of messages over the ROS system.
 */
using namespace std;

class MPCC
{

public:
   MPCC(ros::NodeHandle& nh);
   void run();

private:
   ros::NodeHandle nh_;
   ros::Subscriber odom_sub_;
   string pose_topic;
   void odometry_callback(const nav_msgs::Odometry::ConstPtr &odom_msg);
   void getParameters(ros::NodeHandle& nh);

};


MPCC::MPCC(ros::NodeHandle& nh): nh_(nh)
{
    getParameters(nh_);

    /**
     * The subscribe() call is how you tell ROS that you want to receive messages
     * on a given topic.  This invokes a call to the ROS
     * master node, which keeps a registry of who is publishing and who
     * is subscribing.  Messages are passed to a callback function, here
     * called chatterCallback.  subscribe() returns a Subscriber object that you
     * must hold on to until you want to unsubscribe.  When all copies of the Subscriber
     * object go out of scope, this callback will automatically be unsubscribed from
     * this topic.
     *
     * The second parameter to the subscribe() function is the size of the message
     * queue.  If messages are arriving faster than they are being processed, this
     * is the number of messages that will be buffered up before beginning to throw
     * away the oldest ones.
     */
  odom_sub_ = nh_.subscribe(this->pose_topic, 10, &MPCC::odometry_callback, this);

}

void MPCC::run()
{
  ROS_INFO("Yallah Let go\n");
}

void MPCC::getParameters(ros::NodeHandle& nh)
{
  nh.getParam("pose_topic", pose_topic);
}

void odometry_callback(const nav_msgs::Odometry::ConstPtr &odom_msg)
{
  ROS_INFO("Position x: [%d], y: [%d] ", odom_msg->pose.pose.position.x, odom_msg->pose.pose.position.y);
}


int main(int argc, char **argv){
    ros::init(argc, argv, "mpcc");
    ros::NodeHandle nh;
    MPCC mpcc(nh);
    ros::Rate rate(20);
    while(ros::ok()){
        ros::spinOnce();
        mpcc.run();
        rate.sleep();
    }
    return 0;
}
