#include "placePeopleInZenith.h"

std::vector<double> getVectorFromFile(std::string fname) {
    std::ifstream ifs(fname.c_str());
    std::istream_iterator<double> start(ifs), end;
    std::vector<double> numbers(start, end);
    ifs.close();
    return std::move(numbers);
}

void putVectorInAFile(std::vector<double> &outv, std::string fname) {
    std::ofstream ofs(fname.c_str());
    std::ostream_iterator<double> out_it (ofs, "\n");
    std::copy ( outv.begin(), outv.end(), out_it);
    ofs.close();
}
  
struct Point 
{ 
    double x; 
    double y; 
}; 
  
// Given three colinear points p, q, r, the function checks if 
// point q lies on line segment 'pr' 
bool onSegment(Point p, Point q, Point r) 
{ 
    if (q.x <= std::max(p.x, r.x) && q.x >= std::min(p.x, r.x) && 
            q.y <= std::max(p.y, r.y) && q.y >= std::min(p.y, r.y)) 
        return true; 
    return false; 
} 
  
// To find orientation of ordered triplet (p, q, r). 
// The function returns following values 
// 0 --> p, q and r are colinear 
// 1 --> Clockwise 
// 2 --> Counterclockwise 
int orientation(Point p, Point q, Point r) 
{ 
    double val = (q.y - p.y) * (r.x - q.x) - 
              (q.x - p.x) * (r.y - q.y); 
  
    if (std::abs(val) < epsi) return 0;  // colinear 
    return (val > 0)? 1: 2; // clock or counterclock wise 
} 
  
// The function that returns true if line segment 'p1q1' 
// and 'p2q2' intersect. 
bool doIntersect(Point p1, Point q1, Point p2, Point q2) 
{ 
    // Find the four orientations needed for general and 
    // special cases 
    int o1 = orientation(p1, q1, p2); 
    int o2 = orientation(p1, q1, q2); 
    int o3 = orientation(p2, q2, p1); 
    int o4 = orientation(p2, q2, q1); 
  
    // General case 
    if (o1 != o2 && o3 != o4) 
        return true; 
  
    // Special Cases 
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1 y = (ysup -ymin) / (r4 * std::sin(angle) - ymin) * (y-ymin) + ymin;
    if (o1 == 0 && onSegment(p1, p2, q1)) return true; 
  
    // p1, q1 and p2 are colinear and q2 lies on segment p1q1 
    if (o2 == 0 && onSegment(p1, q2, q1)) return true; 
  
    // p2, q2 and p1 are colinear and p1 lies on segment p2q2 
    if (o3 == 0 && onSegment(p2, p1, q2)) return true; 
  
     // p2, q2 and q1 are colinear and q1 lies on segment p2q2 
    if (o4 == 0 && onSegment(p2, q1, q2)) return true; 
  
    return false; // Doesn't fall in any of the above cases 
} 
  
// Returns true if the point p lies inside the polygon[] with n vertices 
bool isInside(Point polygon[], int n, Point p) 
{ 
    // There must be at least 3 vertices in polygon[] 
    if (n < 3)  return false; 
  
    // Create a point for line segment from p to infinite 
    Point extreme = {INF, p.y}; 
  
    // Count intersections of the above line with sides of polygon 
    int count = 0, i = 0; 
    do
    { 
        int next = (i+1)%n; 
  
        // Check if the line segment from 'p' to 'extreme' intersects 
        // with the line segment from 'polygon[i]' to 'polygon[next]' 
        if (doIntersect(polygon[i], polygon[next], p, extreme)) 
        { 
            // If the point 'p' is colinear with line segment 'i-next', 
            // then check if it lies on segment. If it lies, return true, 
            // otherwise false 
            if (orientation(polygon[i], p, polygon[next]) == 0) 
               return onSegment(polygon[i], p, polygon[next]); 
  
            count++; 
        } 
        i = next; 
    } while (i != 0); 
  
    // Return true if count is odd, false otherwise 
    return count&1;  // Same as (count%2 == 1) 
} 

bool insideDoors(double x, double y, std::vector<double> &doors)
{
    Point p{x,y};
    size_t siz = doors.size();
    size_t nbRectangles = siz / 8;
    for(int i=0; i<nbRectangles; i++){
        Point door[4];
        door[0] = Point{doors[8*i], doors[8*i+1]};
        door[1] = Point{doors[8*i+2], doors[8*i+3]};
        door[2] = Point{doors[8*i+4], doors[8*i+5]};
        door[3] = Point{doors[8*i+6], doors[8*i+7]};
        if (isInside(door, 4, p)) return true;
    }
    return false;
}

void placePeople(std::vector<double> &initX, std::vector<double> &initY)
{
    std::vector<double> doors { getVectorFromFile("../inputs/zenith_up_doors.txt")};
    double r1 = 19.5366802215576;
    double r2 = 34.8401107788086;
    double r3 = 38.601261138916;
    double r4 = 52.2790603637695;
    double maxRandomPosition = 2500;
    int nbRows = 20;
    int nbRowsUp = 15;
    int maxPeoplePerRow = 100;
    double angleMin = 20 * M_PI / 180.0;
    double angleMax = M_PI - angleMin;
    double dangle = (angleMax - angleMin) / maxPeoplePerRow;
    initX.reserve((nbRows + nbRowsUp) * maxPeoplePerRow + maxRandomPosition);
    initY.reserve((nbRows + nbRowsUp) * maxPeoplePerRow + maxRandomPosition);
    for(double radius = r1; radius <= r2; radius += (r2-r1)/nbRows) {
        for(double angle = angleMin; angle <= angleMax; angle += dangle) {
            double x = radius * std::cos(angle);
            double y = radius * std::sin(angle);
            if (!insideDoors(x, y, doors)) {
                initX.push_back(x);
                initY.push_back(y);
            }
        }
    }
    // random position
    
    std::uniform_real_distribution<double> distributionPX(0.0, r1-dr);
    std::uniform_real_distribution<double> distributionTheta(0.0, M_PI);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generatorR (seed);
    std::default_random_engine generatorTheta(seed + 10);
    for(int i=0; i<maxRandomPosition; i++){
        double r = distributionPX(generatorR);
        double theta = distributionTheta(generatorTheta);
        initX.push_back(r * std::cos(theta));
        initY.push_back(r * std::sin(theta));
    }

    // no people on stage

    /* std::uniform_real_distribution<double> distributionSceneX(-9.0, 9.0);
    std::uniform_real_distribution<double> distributionSceneY(-13.0, -1.0);
    for(int i=0; i<10; i++) {
        initX.push_back(distributionSceneX(generatorR));
        initY.push_back(distributionSceneY(generatorTheta));
    } */

    std::vector<double> tmp{getVectorFromFile("../inputs/ZenithWallTop.txt")};
    std::vector<double> xw, yw;
    xw.reserve(tmp.size()/2);
    yw.reserve(tmp.size()/2);
    for(int i=0; i<tmp.size(); i+=2) {
        xw.push_back(tmp[i]);
        yw.push_back(tmp[i+1]);
    }
    std::reverse(xw.begin(), xw.end());
    std::reverse(yw.begin(), yw.end());

    std::vector<double> tmpm{getVectorFromFile("../inputs/zenithTopRmin.txt")};
    std::vector<double> xwm, ywm;
    xwm.reserve(tmpm.size()/2);
    ywm.reserve(tmpm.size()/2);
    for(int i=0; i<tmpm.size(); i+=2) {
        xwm.push_back(tmpm[i]);
        ywm.push_back(tmpm[i+1]);
    }
    std::vector<double> xMax(xw.begin(), xw.begin() + 7);
    std::vector<double> yMax(yw.begin(), yw.begin() + 7);
    std::vector<double> xMin(xwm.begin() + 7, xwm.end());
    std::vector<double> yMin(ywm.begin() + 7, ywm.end());
    std::reverse(xMin.begin(), xMin.end());
    std::reverse(yMin.begin(), yMin.end());

    //boost::math::barycentric_rational<double> s(xw.data(), yw.data(), xw.size());
    //boost::math::barycentric_rational<double> sm(xwm.data(), ywm.data(), xwm.size());
    linearInterp s(xw, yw);
    linearInterp sm(xwm, ywm);
    linearInterp sx(yMax, xMax);
    linearInterp sxm (yMin, xMin);
    angleMax = M_PI * 0.5 - dangle;
    for(double radius = r3; radius <= r4 - dr; radius += (r4-dr-r3)/nbRowsUp) {
        for(double angle = angleMin; angle <= angleMax; angle += dangle) {
            double x = radius * std::cos(angle);
            double y = radius * std::sin(angle);
            double ysup = s(x);
            double ymin = sm(x);
            double xsup = sx(y);
            double xmin = sxm(y);
            y = (ysup -ymin) / ((r4 - r3) * std::sin(angle)) * (y - r3 * std::sin(angle)) + ymin;
            //x = (xsup - xmin) / (r4 * std::cos(angle) - r3 * std::cos(angle)) * (x - r3 * std::cos(angle)) + xmin;
            //x = xsup;
            //y = ysup;
            if (!insideDoors(x, y, doors) && x < xw.back()) {
                initX.push_back(x);
                initX.push_back(-x);
                initY.push_back(y);
                initY.push_back(y);
            }
        }
    }

    //putVectorInAFile(initX, "initX_zenith_middle.dat");
    //putVectorInAFile(initY, "initY_zenith_middle.dat");
}
