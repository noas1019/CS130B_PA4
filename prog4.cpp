#include <iostream>
#include <vector>
#include <sstream>
#include <random>
#include <set>
#include <math.h>
#include <utility>
#include <float.h>
#include <algorithm>
#include <cmath>
using namespace std;

// custom Vertex class to hold x,y coordinates
class Vertex {
public:
    double x;
    double y;
    Vertex();
    Vertex(double x_cord, double y_cord);
};

// default constructor
Vertex::Vertex() {
    this->x = 0.0;
    this->y = 0.0;
}

// parameterized constructor
Vertex::Vertex(double x_cord, double y_cord) {
    this->x = x_cord;
    this->y = y_cord;
}

// simple function to find median value in vector of doubles
double findMedian(vector<double> &errors) {
    int index = errors.size() / 2;
    nth_element(errors.begin(), errors.begin() + index, errors.end());
    return errors[index];
}

// function will calculate shortest distance from Vertex A to best-fit line formed by Vertices Start and End
double minDistance(Vertex A, Vertex Start, Vertex End) {
    
    // difference between start and end vertices
    Vertex Start_End;
    Start_End.x = End.x - Start.x;
    Start_End.y = End.y - Start.y;
    
    // difference between end and target vertices
    Vertex End_Target;
    End_Target.x = A.x - End.x;
    End_Target.y = A.y - End.y;
    
    // difference between start and target vertices
    Vertex Start_Target;
    Start_Target.x = A.x - Start.x;
    Start_Target.y = A.y - Start.y;

    // calculate dot products of Start_End, End_Target and Start_End, Start_Target
    double Start_End_End_Target = (Start_End.x * End_Target.x + Start_End.y * End_Target.y);
    double Start_End_Start_Target = (Start_End.x * Start_Target.x + Start_End.y * Start_Target.y);

    // calculate distance for different cases
    double result = 0;
    if (Start_End_End_Target > 0) {
        double xx = A.x - End.x;
        double yy = A.y - End.y;
        result = sqrt(xx * xx + yy * yy);
    }
    else if (Start_End_Start_Target < 0) {
        double xx = A.x - Start.x;
        double yy = A.y - Start.y;
        result = sqrt(xx * xx+ yy * yy);
    }
    else {
        double x1 = Start_End.x;
        double y1 = Start_End.y;
        double x2 = Start_Target.x;
        double y2 = Start_Target.y;
        double divisor = sqrt(x1 * x1 + y1 * y1);
        result = abs(x1 * y2 - y1 * x2) / divisor;
    }
    return result;
}

// linear least squares method to minimize error and return best-fit a and b coefficients
pair<double, double> LeastSquares(vector<Vertex> &points, double min_y) {
    // initialize variables
    double top, bottom, xAVG, yAVG, A, B;
    double vectorSize= points.size();

    // if vector is singular we just return y-intercept
    if (vectorSize == 1) {
        A = 0.0;
        B = min_y;
        pair<double, double> coefficients;
        coefficients.first = A;
        coefficients.second = B;
        return coefficients;
    }
    
    // calculate x and y averages
    xAVG = 0.0;
    yAVG = 0.0;
    for (int i = 0; i < vectorSize; i++) {
        xAVG = xAVG + points[i].x;
        yAVG = yAVG + points[i].y;
    }
    xAVG = xAVG / vectorSize;
    yAVG = yAVG / vectorSize;
    
    // calculate top and bot
    top = 0.0;
    bottom = 0.0;
    for (int i = 0; i < vectorSize; i++) {
        top = top + (points[i].x - xAVG) * (points[i].y - yAVG);
        bottom = bottom + (points[i].x - xAVG) * (points[i].x - xAVG);
    }

    // calculate resulting coefficients
    A = top / bottom;
    B = yAVG - A * xAVG;
    pair<double, double> coefficients;
    coefficients.first = A;
    coefficients.second = B;
    return coefficients;
}


// driver code
int main() {
    
    // read file input into vector of vertices
	vector<Vertex> points;
    string vertex;
    while (getline(cin, vertex)) {
        double x, y;
        istringstream iss(vertex);
        iss >> x;
        iss >> y;
        Vertex point(x, y);
        points.push_back(point);
    }

    // we are doing 2*log(n) iterations
    const int iterations = 2*log(points.size());
    int rands1[iterations] = {};
    int rands2[iterations] = {};

    // this is a C++ random engine that should give a pretty random integer within the range of our vector
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0, points.size()-1);

    // initialize our solution variables
    Vertex best1;
    Vertex best2;
    double a, b;

    // set min_median to highest double value initially
    double min_median = DBL_MAX;
    
    // Iterating 2*log(n) times
    for (int i = 0; i < iterations; ++i) {

        // Use uniform distribution generator to get random vertices
        rands1[i] = distribution(generator);
        rands2[i] = distribution(generator);
        Vertex random1 = points[rands1[i]];
        Vertex random2 = points[rands2[i]];

        // Calculate slope and y-intercept with random vertices
        double slope = (random2.y - random1.y) / (random2.x-random1.x);
        double y_intercept = random1.y - (slope * random1.x);
        
        // Store vector of distances for each point from best fit line
        vector<double> errors;
        for (Vertex v : points) {
            double error = minDistance(v, random1, random2);
            errors.push_back(error);
        }
        
        // Find median error in vector of errors
        double median = findMedian(errors);
        if (median < min_median) {
            min_median = median;
            a = slope;
            b = y_intercept;
            best1 = random1;
            best2 = random2;
        }
    }

    // loop through original vertex of points and calculate minDistance to best-fit line to see if they fall under min_median, if so, add to another vertex of vertices we will use as samples for least-square regression
    // we also calculate minimum y value among fitting points for edge cases in LeastSquares
    vector<Vertex> fitPoints;
    double min_Y = DBL_MAX;
    for (Vertex v : points) {
        double vError = minDistance(v, best1, best2);
        if (vError < min_median) {
            if (v.y < min_Y) {
                min_Y = v.y;
            }
            fitPoints.push_back(v);
        }
    }
    
    // output least-squares best-fit coefficients to console rounded to 2 decimal places
    pair<double, double> result = LeastSquares(fitPoints, min_Y);
    double A = result.first;
    double B = result.second;
    cout << ceil(A * 100.0) / 100.0 << " " << ceil(B * 100.0) / 100.0 << endl;
    
    return 0;
}