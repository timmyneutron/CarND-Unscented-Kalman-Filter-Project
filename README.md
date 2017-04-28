![Udacity - Self-Driving Car NanoDegree](https://s3.amazonaws.com/udacity-sdc/github/shield-carnd.svg)

# Udacity Self-Driving Car Nanodegree: Term 2
# Project #2: Unscented Kalman Filter

## Introduction
This is a project for Udacity's Self-Driving Car Nanodegree. It implements an Unscented Kalman Filter to estimate a car's state based on radar and laser measurements.


## Concepts and Classes
Concepts explored in this project:

  - Unscented Kalman Filters
  - Radar and laser sensor fusion
  - Non-linear process models like CTRV (Constant Turn Rate and Velocity)

Relevant classes:

  - [Artificial Intelligence for Robotics](https://classroom.udacity.com/courses/cs373)

## Getting Started
Code is stored in the `src` directory.

Input data and output data (state estimations and final root mean squared error) are in the `data` directory.

To run the code, navigate to the `build` directory and execute the following:

`./UnscentedKF ../data/obj_pose-laser-radar-synthetic-input.txt ../data/output.txt`

or

`./UnscentedKF ../data/obj_pose-laser-radar-synthetic-input.txt ../data/output.txt`
