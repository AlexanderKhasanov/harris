#include <fstream>
#include <string>
#include <boost/filesystem.hpp>
#include "detector_harrisa.hpp"


using namespace std;
namespace fs = boost::filesystem;


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

double comparison_descriptor(const vector<descriptor_type>& img1, const vector<descriptor_type>& img2)
{
    unsigned int number_similar = 0;
    unsigned int number_varioous = 0;
    unsigned int max_distance = 105;

    for( auto dis1 = img1.begin(); dis1 != img1.end(); ++dis1 )
    {
        bool q = false;
        for( auto dis2 = img2.begin(); dis2 != img2.end(); ++dis2 )
        {

            int k = hamming_distance(*dis1, *dis2);
            if ( hamming_distance(*dis1, *dis2) <= max_distance)
            {
                ++number_similar;
                q = true;
                break;
            }
        }
        if (!q)
            ++number_varioous;
    }
    if ( number_varioous == 0 )
        number_varioous = 1;
    return double(number_similar)/number_varioous;
}

int main(int argc, char *argv[])
{
    ifstream file(argv[1]);
    //ifstream file("/home/alexander/study/Алгоритмы/Harris/input.txt");


    if(!file.is_open())
    {
        return 0;
    }

    detector_harrisa detecter;

    vector<string> result;
    vector<string> command;
    string tmp;

    while ( getline(file, tmp) )
    {
        if ( tmp != "" )
            command.push_back(tmp);
    }

    if ( command.size() != 3 )
        return 0;

    string pattern_path = command[0];
    string images_path = command[1];
    string result_path = command[2];

    matrix_type pattern_img_pixel = imread( pattern_path, 0 );
    Mat_<Vec3b> img = imread(pattern_path, 1);

    vector<key_point> key_point1 = detecter.search_key_point(pattern_img_pixel);

    vector<descriptor_type> pattern_descriptor = detecter.search_descriptor( pattern_img_pixel );

    fs::directory_iterator end_itr;

    for ( fs::directory_iterator it( images_path ); it != end_itr; ++it )
    {
        string img_path = it->path().string();

        matrix_type img_pixel = imread(img_path, 0);

        vector<descriptor_type> img_descriptor = detecter.search_descriptor( img_pixel );

        double k = comparison_descriptor( pattern_descriptor, img_descriptor );

        if ( k >= 1)
        {
            result.push_back( img_path );
            //cout<<" (Уверенность - "<<k<<")"<<endl;
        }
    }

    ofstream file_result( result_path );

    if (!file_result.is_open())
        return 0;

    for (const auto & it : result)
    {
        file_result << it << endl;
    }
    file_result.close();

//    for(auto current = key_point1.begin(); current != key_point1.end(); ++current)
//    {
//        Point a(current->coordinates.x, current->coordinates.y);
//        circle(img, a, 3, Scalar(0, 0, 255), 1);
//    }
//
//    imwrite("/home/alexander/study/Алгоритмы/Harris/resultK11111.jpeg", img);

    return 0;
}