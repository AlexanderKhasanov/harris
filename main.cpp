#include "detector_harrisa.hpp"

int hamming_distance(const vector<bool>& vec1, const vector<bool>& vec2)
{
    if( vec1.size() != vec2.size() )
        return -1;
    int distance = 0;
    for(size_type i = 0; i < vec1.size(); ++i)
    {
        if( vec1[i] != vec2[i] )
            ++distance;
    }
    return distance;
}


int main() {
    // получаем картинку
    Mat_<double> image = imread( "/home/alexander/study/Алгоритмы/Harris/1.jpg", 0);

    Mat_<Vec3b> image2 = imread( "/home/alexander/study/Алгоритмы/Harris/1.jpg", 1);

    detector_harrisa detecter1;

    vector<key_point> key_point1 = detecter1.search_key_point(image);

    vector<descriptor_type> descriptor1 = detecter1.brief(image, key_point1);

    auto it = descriptor1.begin();

    vector<bool> vec = it->second;

    for (auto it = vec.begin(); it != vec.end(); ++it)
        cout<< *it<< " ";


//    for(auto current = points.begin(); current != points.end(); ++current)
//    {
//        cout<<"x = "<< current->x<<endl;
//        cout<<"y = "<< current->y<<endl;
//        cout<<"intens = "<< current->intensity<<"\n"<<endl;
//        Point a(current->x, current->y);
//        circle(image2, a, 3, Scalar(0, 0, 255), 1);
//    }
//
//    imwrite("/home/alexander/study/Алгоритмы/Harris/result12.jpeg", image2);
//    cout << "s = "<< points.size()<<endl;

    return 0;
}