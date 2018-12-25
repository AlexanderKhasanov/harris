#include <iostream>
#include <cstdlib>
#include <opencv2/opencv.hpp>
#include <cmath>
#include <vector>
#include <ctime>

using namespace cv;
using namespace std;

struct point
{
    int x;
    int y;

    point(int _x, int _y)
        : x (_x )
        , y (_y )
    {}

    ~point() = default;
};

struct key_point
{
    point coordinates;
    double intensity;

    key_point(int x, int y, double intens)
        : coordinates(x, y)
        , intensity(intens)
    {}
};

using value_type = double;
using matrix_type = Mat_<value_type>;
using size_type = size_t;
using descriptor_type = vector<bool>;

class detector_harrisa
{
    //параметры детектора
    double _k;
    size_type _bobl_raduis; //радиус в котором выбирается точка с наибольшей интенсивностью
    int _min_intensiti; //минимальная интенсивность
    size_type _window_size; //размер окна вычисления интенсивности

    //параметры дескриптора
    size_type _dimensionality; //размер ветора сравнений
    size_type _size_area_point; //размер области вокруг особой точки
    size_type _area_for_average; //размер окна для вычисления средней интенсивности для сравнения

public:
    detector_harrisa();

    detector_harrisa(size_type window_size, int min_intensiti, size_type bobl_raduis, double k,
                    size_type dimensionality, size_type size_area_point, size_type area_for_average);

    ~detector_harrisa() = default;

    std::vector<key_point> search_key_point(const matrix_type& image);

    vector<descriptor_type> brief(const matrix_type & pixel_image, const std::vector<key_point>& key_points);

    vector<descriptor_type> search_descriptor(const matrix_type& image);


private:
    void normalize(matrix_type& image);
    matrix_type gradient_x(const matrix_type& image);
    matrix_type gradient_y(const matrix_type& image);

    double average_intensity (int x, int y, const matrix_type & pixel_image) const;

    size_type window_size() const;

    double k() const;

    size_type bobl_raduis() const;

    int min_intensiti() const;

    size_type dimensionality() const;

    size_type size_area_point() const;

    size_type area_for_average() const;
};