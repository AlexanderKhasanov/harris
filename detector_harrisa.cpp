#include "detector_harrisa.hpp"

detector_harrisa::detector_harrisa()
    : _window_size (5)
    , _min_intensiti (10)
    , _bobl_raduis (5)
    , _k (0.05)
    , _dimensionality (256)
    , _size_area_point (30)
    , _area_for_average (5)
{}

detector_harrisa::detector_harrisa(size_type window_size, int min_intensiti, size_type bobl_raduis, double k,
                                   size_type dimensionality, size_type size_area_point, size_type area_for_average)
        : _window_size ( window_size )
        , _min_intensiti ( min_intensiti )
        , _bobl_raduis ( bobl_raduis )
        , _k ( k )
        , _dimensionality ( dimensionality )
        , _size_area_point ( size_area_point )
        , _area_for_average ( area_for_average )
{}

size_type detector_harrisa::window_size() const
{
    return _window_size;
}

double detector_harrisa::k() const
{
    return _k;
}

size_type detector_harrisa::bobl_raduis() const
{
    return _bobl_raduis;
}

int detector_harrisa::min_intensiti() const
{
    return _min_intensiti;
}

size_type detector_harrisa::dimensionality() const
{
    return _dimensionality;
}

size_type detector_harrisa::size_area_point() const
{
    return _size_area_point;
}

size_type detector_harrisa::area_for_average() const
{
    return _area_for_average;
}

void detector_harrisa::normalize(matrix_type& image)
{
    for(size_type h = 0; h < image.rows; ++h)
        for(size_type w = 0; w < image.cols; ++w)
            image(h, w) /= 255;
}

matrix_type detector_harrisa::gradient_x(const matrix_type &image)
{
    size_type kernel_size = 3;
    vector<vector<int>> kernel_x = { {-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1} };

    matrix_type I_x(image.rows, image.cols);
    for(size_type h = 1; h < image.rows - 1; ++h)
        for(size_type w = 1; w < image.cols - 1; ++w)
        {
            double sum = 0;
            for(size_type i = 0; i < kernel_size; ++i)
                for(size_type j = 0; j < kernel_size; ++j)
                {
                    double current_pixel = image( h-1 + i, w-1 + j);
                    sum += current_pixel * kernel_x[i][j];
                }
            I_x( h, w ) = sum;
        }
    return I_x;
}

matrix_type detector_harrisa::gradient_y(const matrix_type &image)
{
    size_type kernel_size = 3;
    vector<vector<int>> kernel_x = { {-1, -2, -1}, {0, 0, 0}, {1, 2, 1} };

    matrix_type I_y(image.rows, image.cols);

    for(size_type h = 1; h < image.rows - 1; ++h)
        for(size_type w = 1; w < image.cols - 1; ++w)
        {
            double sum = 0;
            for(size_type i = 0; i < kernel_size; ++i)
                for(size_type j = 0; j < kernel_size; ++j)
                {
                    double current_pixel = image( h-1 + i, w-1 + j);
                    sum += current_pixel * kernel_x[i][j];
                }
            I_y( h, w ) = sum;
        }
    return I_y;
}


std::vector<key_point> detector_harrisa::search_key_point(const matrix_type &image)
{
    std::vector<key_point> key_points;
    bool add = false;
    size_type offset = window_size()/2;//?

    matrix_type normalize_image;
    image.copyTo(normalize_image);
    normalize(normalize_image);

    Mat_<double> I_x = gradient_x(normalize_image);
    Mat_<double> I_y = gradient_y(normalize_image);

    Mat_<double> I_xx;
    Mat_<double> I_xy;
    Mat_<double> I_yy;

    pow(I_x, 2.0, I_xx);
    pow(I_y, 2.0, I_yy);
    multiply(I_x, I_y, I_xy);

    for (size_t h = 1 + offset; h < normalize_image.rows - offset - 1; ++h)
        for (size_t w = 1 + offset; w < normalize_image.cols - offset - 1; ++w)
        {
            double sum_xx = 0;
            double sum_xy = 0;
            double sum_yy = 0;
            for (size_t i = 0; i < window_size(); ++i)
                for(size_t j = 0; j < window_size(); ++j)
                {
                    sum_xx += I_xx(h - offset + i, w - offset + j);
                    sum_xy += I_xy(h - offset + i, w - offset + j);
                    sum_yy += I_yy(h - offset + i, w - offset + j);
                }
            double det = sum_xx * sum_yy - sum_xy * sum_xy;
            double trace = sum_xx + sum_yy;
            double intensty = det - k() * trace * trace ;

            if ( intensty > min_intensiti() )
            {
                if (key_points.size() == 0)
                {
                    key_points.push_back( key_point( w, h, intensty ) );
                }
                else
                {
                    for(auto current = key_points.begin(); current != key_points.end(); ++current)
                    {
                        int delta_x = w - current->coordinates.x;
                        int delta_y = h - current->coordinates.y;
                        int rad = pow(delta_x, 2.0) + pow( delta_y, 2.0 );
                        if( rad <= pow( bobl_raduis(), 2.0 ) )
                        {
                            add = false;
                            if ( intensty > current->intensity )
                            {
                                current->coordinates.x = w;
                                current->coordinates.y = h;
                                current->intensity = intensty;
                            }
                            break;
                        }
                        else
                            add = true;
                    }
                    if(add)
                        key_points.push_back( key_point( w, h, intensty ) );
                }
            }

        }
    return key_points;
}

double detector_harrisa::average_intensity(int x, int y, const matrix_type & pixel_image) const
{
    size_type offset = area_for_average()/2;
    int sum = 0;
    int count = 0;
    for (int h = 0; h < area_for_average(); ++h)
        for (int w = 0; w < area_for_average(); ++w)
        {
            sum += pixel_image( y - offset + h, x - offset + w);
            ++count;
        }
    return double(sum)/count;
}

vector<descriptor_type> detector_harrisa::brief(const matrix_type &pixel_image,
                                                const std::vector<key_point> &key_points)
{
    srand(time(0));
    vector<descriptor_type> descriptors;
    for (auto current_key = key_points.begin(); current_key != key_points.end(); ++ current_key)
    {
        vector<bool> vec;
        for (size_type i = 0; i < _dimensionality; ++i)
        {
            size_type offset = size_area_point() / 2;
            int q1 = current_key->coordinates.x + offset;
            int q2 = current_key->coordinates.x - offset;
            int q3 = current_key->coordinates.y + offset;
            int q4 = current_key->coordinates.y - offset;
            while ( q1 >= pixel_image.cols || q2 < 0 || q3 >= pixel_image.rows || q4 < 0 )
            {
                offset /= 2;
                q1 = current_key->coordinates.x + offset;
                q2 = current_key->coordinates.x - offset;
                q3 = current_key->coordinates.y + offset;
                q4 = current_key->coordinates.y - offset;
            }

            int delta_x1 = 0;
            int delta_x2 = 0;

            int delta_y1 = 0;
            int delta_y2 = 0;

            size_type interval = offset - area_for_average()/2;
            if ( interval > 0 )
            {
                delta_x1 = -interval + rand() % ( 2 * interval );
                delta_x2 = -interval + rand() % ( 2 * interval );

                delta_y1 = -interval + rand() % ( 2 * interval );
                delta_y2 = -interval + rand() % ( 2 * interval );
            }


            double I1 = average_intensity(current_key->coordinates.x + delta_x1,
                                          current_key->coordinates.y + delta_y1, pixel_image);

            double I2 = average_intensity(current_key->coordinates.x + delta_x2,
                                          current_key->coordinates.y + delta_y2, pixel_image);

            if (I1 < I2)
                vec.push_back(1);
            else
                vec.push_back(0);
        }
        descriptors.push_back( vec );
    }
    return descriptors;
}

vector<descriptor_type> detector_harrisa::search_descriptor(const matrix_type &image)
{
    return brief( image, search_key_point(image) );
}