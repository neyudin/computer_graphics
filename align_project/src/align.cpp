#include "align.h"
#include <string>
#include "io.h"//for debugging
#include <cmath>
#include <typeinfo>

using std::string;
using std::cout;
using std::endl;

using std::tuple;
using std::tie;
using std::make_tuple;

struct Pos
{
    int row_shift, col_shift, rows, cols;
    Pos(int r = 0, int c = 0, int rs = 0, int cs = 0): row_shift(r), col_shift(c), rows(rs), 
    cols(cs) {}
};

Image find_mse(Image &srcA, Image &srcB, int delta_shift, uint channelA)
{
    Pos find_pos;
    double min_error = 0.0;
    for (int i = -delta_shift; i < delta_shift + 1; ++i) {
        int rows = std::min(srcA.n_rows, srcB.n_rows);
        int dr = std::abs(static_cast<int>(srcA.n_rows) - static_cast<int>(srcB.n_rows));
        if (srcA.n_rows >= srcB.n_rows) {
            if (i < 0) {
                rows += i;
            }
            if (i > dr) {
                rows += dr - i;
            }
        } else {
            if (i > 0) {
                rows -= i;
            }
            if (i < -dr) {
                rows += dr + i;
            }
        }
        uint a_row_shift = (i > 0) ? i : 0;
        uint b_row_shift = ((i < 0) ? i : 0);
        for (int j = -delta_shift; j < delta_shift + 1; ++j) {
            double error = 0.0;
            unsigned long long sum_error = 0;
            int cols = std::min(srcA.n_cols, srcB.n_cols);
            int dc = std::abs(static_cast<int>(srcA.n_cols) - static_cast<int>(srcB.n_cols));
            if (srcA.n_cols >= srcB.n_cols) {
                if (j < 0) {
                    cols += j;
                }
                if (j > dc) {
                    cols += dc - j;
                }
            } else {
                if (j > 0) {
                    cols -= j;
                }
                if (j < -dc) {
                    cols += dc + j;
                }
            }
            uint a_col_shift = (j > 0) ? j : 0;
            uint b_col_shift = (j < 0) ? j : 0;
            for (int row = 0; row < rows; ++row) {
                for (int col = 0; col < cols; ++col) {
                    uint a_candidate = std::get<1>(srcA(row + a_row_shift, col + a_col_shift));
                    uint b_candidate = std::get<1>(srcB(row - b_row_shift, col - b_col_shift));
                    sum_error += (a_candidate - b_candidate) * (a_candidate - b_candidate);                    
                }
            }
            error = static_cast<double>(sum_error) / (rows * cols);
            if ((i == -delta_shift) && (j == -delta_shift)) {
                min_error = error;
                find_pos.row_shift = i;
                find_pos.col_shift = j;
                find_pos.rows = rows;
                find_pos.cols = cols;
            } else {
                if (error < min_error) {
                    min_error = error;
                    find_pos.row_shift = i;
                    find_pos.col_shift = j;
                    find_pos.rows = rows;
                    find_pos.cols = cols;
                }
            }
        }
    }
    Image result(find_pos.rows, find_pos.cols);
    int a_row_shift = ((find_pos.row_shift > 0) ? find_pos.row_shift : 0);
    int b_row_shift = ((find_pos.row_shift < 0) ? find_pos.row_shift : 0);
    int a_col_shift = ((find_pos.col_shift > 0) ? find_pos.col_shift : 0);
    int b_col_shift = ((find_pos.col_shift < 0) ? find_pos.col_shift : 0);
    switch (channelA) {
        case 2: {//srcA = R, srcB = G
            for (int row = 0; row < find_pos.rows; ++row) {
                for (int col = 0; col < find_pos.cols; ++col) {
                    std::get<0>(result(row, col)) = std::get<0>(srcA(row + a_row_shift, 
                        col + a_col_shift));
                    std::get<1>(result(row, col)) = std::get<1>(srcB(row - b_row_shift, 
                        col - b_col_shift));
                }
            }
            break;
        }
        case 0: {//srcA = B, srcB = RG
            for (int row = 0; row < find_pos.rows; ++row) {
                for (int col = 0; col < find_pos.cols; ++col) {
                    std::get<0>(result(row, col)) = std::get<0>(srcB(row - b_row_shift, 
                        col - b_col_shift));
                    std::get<1>(result(row, col)) = std::get<1>(srcB(row - b_row_shift, 
                        col - b_col_shift));
                    std::get<2>(result(row, col)) = std::get<2>(srcA(row + a_row_shift, 
                        col + a_col_shift));
                }
            }
            break;
        }
        default: break;
    }
    return result;
}

struct RGB
{
    Image R, G, B;
    RGB(Image r, Image g, Image b): R(r), G(g), B(b) {}
};

struct FindMax
{
    int hist_max, max1, max2;
    FindMax(uint hm = 0, uint m1 = 0, uint m2 = 0): hist_max(hm), max1(m1), max2(m2) {}
};

struct Edges
{
    int rows[6], cols[2];
    Edges(int r[6] = {0}, int c[2] = {0}) {
        for (int i = 0; i < 6; ++i) {
            rows[i] = r[i];
        }
        for (int i = 0; i < 2; ++i) {
            cols[i] = c[i];
        }
    }
};

RGB cut_image(Image &srcImage, Image &edge_map, double cut_scale) {
    FindMax find_max_side1, find_max_side2;
    Edges edges;
    int eps_h = cut_scale * std::min(srcImage.n_cols, srcImage.n_rows);
    int eps_w = cut_scale * std::min(srcImage.n_cols, srcImage.n_rows);
    //horizontal cut
    for (int j = 0; j < 2; ++j) {
        for (int i = 0; i < eps_h + 1; ++i) {
            int sum_up = 0, sum_down = 0;
            for (int col = 0; col < static_cast<int>(edge_map.n_cols); ++col) {
                sum_up += std::get<0>(edge_map(i, col));
                sum_down += std::get<0>(edge_map(edge_map.n_rows - i - 1, col));
            }
            if (sum_up > find_max_side1.hist_max) {
                find_max_side1.hist_max = sum_up;
                if (j == 0) {
                    find_max_side1.max1 = i;
                } else {
                    find_max_side1.max2 = i;
                }
            }
            if (sum_down > find_max_side2.hist_max) {
                find_max_side2.hist_max = sum_down;
                if (j == 0) {
                    find_max_side2.max1 = edge_map.n_rows - i - 1;
                } else {
                    find_max_side2.max2 = edge_map.n_rows - i - 1;
                }
            }
        }
        if (j == 0) {
            find_max_side1.hist_max = 0;
            find_max_side2.hist_max = 0;
            if (find_max_side1.max1 > 0) {
                for (int col = 0; col < static_cast<int>(edge_map.n_cols); ++col) {
                    std::get<0>(edge_map(find_max_side1.max1 - 1, col)) = 0;
                }
                if (find_max_side1.max1 > 1) {
                    for (int col = 0; col < static_cast<int>(edge_map.n_cols); ++col) {
                        std::get<0>(edge_map(find_max_side1.max1 - 2, col)) = 0;
                    }
                }
            }
            if (find_max_side2.max1 < static_cast<int>(edge_map.n_rows) - 1) {
                for (int col = 0; col < static_cast<int>(edge_map.n_cols); ++col) {
                    std::get<0>(edge_map(find_max_side2.max1 + 1, col)) = 0;
                }
                if (find_max_side2.max1 < static_cast<int>(edge_map.n_rows) - 2) {
                    for (int col = 0; col < static_cast<int>(edge_map.n_cols); ++col) {
                        std::get<0>(edge_map(find_max_side2.max1 + 2, col)) = 0;
                    }
                }
            }
            for (int shift = 0; shift < 3; ++shift) {
                for (int col = 0; col < static_cast<int>(edge_map.n_cols); ++col) {
                    std::get<0>(edge_map(find_max_side1.max1 + shift, col)) = 0;
                    std::get<0>(edge_map(find_max_side2.max1 - shift, col)) = 0;
                }
            }
        }
    }
    edges.rows[0] = std::max(find_max_side1.max1, find_max_side1.max2);
    edges.rows[5] = std::min(find_max_side2.max1, find_max_side2.max2);
    //vertical cut
    find_max_side1 = {0, 0, 0}, find_max_side2 = {0, 0, 0};
    for (int j = 0; j < 2; ++j) {
        for (int i = 0; i < eps_w + 1; ++i) {
            int sum_left = 0, sum_right = 0;
            for (int row = 0; row < static_cast<int>(edge_map.n_rows); ++row) {
                sum_left += std::get<0>(edge_map(row, i));
                sum_right += std::get<0>(edge_map(row, edge_map.n_cols - i - 1));
            }
            if (sum_left > find_max_side1.hist_max) {
                find_max_side1.hist_max = sum_left;
                if (j == 0) {
                    find_max_side1.max1 = i;
                } else {
                    find_max_side1.max2 = i;
                }
            }
            if (sum_right > find_max_side2.hist_max) {
                find_max_side2.hist_max = sum_right;
                if (j == 0) {
                    find_max_side2.max1 = edge_map.n_cols - i - 1;
                } else {
                    find_max_side2.max2 = edge_map.n_cols - i - 1;
                }
            }
        }
        if (j == 0) {
            find_max_side1.hist_max = 0;
            find_max_side2.hist_max = 0;
            if (find_max_side1.max1 > 0) {
                for (int row = 0; row < static_cast<int>(edge_map.n_rows); ++row) {
                    std::get<0>(edge_map(row, find_max_side1.max1 - 1)) = 0;
                }
                if (find_max_side1.max1 > 1) {
                    for (int row = 0; row < static_cast<int>(edge_map.n_rows); ++row) {
                        std::get<0>(edge_map(row, find_max_side1.max1 - 2)) = 0;
                    }
                }
            }
            if (find_max_side2.max1 < static_cast<int>(edge_map.n_cols) - 1) {
                for (int row = 0; row < static_cast<int>(edge_map.n_rows); ++row) {
                    std::get<0>(edge_map(row, find_max_side2.max1 + 1)) = 0;
                }
                if (find_max_side2.max1 < static_cast<int>(edge_map.n_cols) - 2) {
                    for (int row = 0; row < static_cast<int>(edge_map.n_rows); ++row) {
                        std::get<0>(edge_map(row, find_max_side2.max1 + 2)) = 0;
                    }
                }
            }
            for (int shift = 0; shift < 3; ++shift) {
                for (int row = 0; row < static_cast<int>(edge_map.n_rows); ++row) {
                    std::get<0>(edge_map(row, find_max_side1.max1 + shift)) = 0;
                    std::get<0>(edge_map(row, find_max_side2.max1 - shift)) = 0;
                }
            }
        }
    }
    edges.cols[0] = std::min(std::max(find_max_side1.max1, find_max_side1.max2), 2 * eps_w / 3);
    edges.cols[1] = std::max(std::min(find_max_side2.max1, find_max_side2.max2), 
        static_cast<int>(srcImage.n_cols) - 2 * eps_w / 3);
    int res_height = edges.rows[5] - edges.rows[0];
    int pred_row1 = edges.rows[0] + res_height / 3;
    int pred_row2 = edges.rows[0] + 2 * (res_height / 3);
    //cut between B and G and R
    find_max_side1 = {0, 0, 0}, find_max_side2 = {0, 0, 0};
    for (int j = 0; j < 2; ++j) {
        for (int i = -eps_h; i < eps_h + 1; ++i) {
            int sum_BG = 0, sum_GR = 0;
            for (int col = 0; col < static_cast<int>(edge_map.n_cols); ++col) {
                sum_BG += std::get<0>(edge_map(pred_row1 + i, col));
                sum_GR += std::get<0>(edge_map(pred_row2 + i, col));
            }
            if (sum_BG > find_max_side1.hist_max) {
                find_max_side1.hist_max = sum_BG;
                if (j == 0) {
                    find_max_side1.max1 = pred_row1 + i;
                } else {
                    find_max_side1.max2 = pred_row1 + i;
                }
            }
            if (sum_GR > find_max_side2.hist_max) {
                find_max_side2.hist_max = sum_GR;
                if (j == 0) {
                    find_max_side2.max1 = pred_row2 + i;
                } else {
                    find_max_side2.max2 = pred_row2 + i;
                }
            }
        }
        if (j == 0) {
            find_max_side1.hist_max = 0;
            find_max_side2.hist_max = 0;
            for (int shift = -2; shift < 3; ++shift) {
                for (int col = 0; col < static_cast<int>(edge_map.n_cols); ++col) {
                    std::get<0>(edge_map(find_max_side1.max1 + shift, col)) = 0;
                    std::get<0>(edge_map(find_max_side2.max1 + shift, col)) = 0;
                }
            }
        }
    }
    edges.rows[1] = std::min(find_max_side1.max1, find_max_side1.max2);
    edges.rows[2] = std::max(find_max_side1.max1, find_max_side1.max2);
    edges.rows[3] = std::min(find_max_side2.max1, find_max_side2.max2);
    edges.rows[4] = std::max(find_max_side2.max1, find_max_side2.max2);
    Image B, G, R;
    if ((std::abs(edges.rows[1] - edges.rows[2]) < 2 * eps_h / 3) && 
        (std::abs(edges.rows[3] - edges.rows[4]) < 2 * eps_h / 3) && (edges.rows[1] != 0) && 
        (edges.rows[2] != 0) && (edges.rows[3] != 0) && (edges.rows[4] != 0)) {
        B = srcImage.submatrix(edges.rows[0], edges.cols[0], edges.rows[1] - edges.rows[0], 
            edges.cols[1] - edges.cols[0]);
        G = srcImage.submatrix(edges.rows[2], edges.cols[0], edges.rows[3] - edges.rows[2], 
            edges.cols[1] - edges.cols[0]);
        R = srcImage.submatrix(edges.rows[4], edges.cols[0], edges.rows[5] - edges.rows[4], 
            edges.cols[1] - edges.cols[0]);
    } else {
        B = srcImage.submatrix(0, edges.cols[0], srcImage.n_rows / 3, 
            edges.cols[1] - edges.cols[0]);
        G = srcImage.submatrix(srcImage.n_rows / 3, edges.cols[0], srcImage.n_rows / 3, 
            edges.cols[1] - edges.cols[0]); 
        R = srcImage.submatrix(2 * (srcImage.n_rows / 3), edges.cols[0], 
            srcImage.n_rows - 2 * (srcImage.n_rows / 3), edges.cols[1] - edges.cols[0]);
    }
    return RGB(R, G, B);
}



Image mirroring(Image &src_image, uint radius) {
    Image mirrored(src_image.n_rows + 2 * radius, src_image.n_cols + 2 * radius);
    uint end_row = mirrored.n_rows - radius;
    uint end_col = mirrored.n_cols - radius;
    for (uint row = radius; row < end_row; ++row) {
        for (uint col = radius; col < end_col; ++col) {
            std::get<0>(mirrored(row, col)) = std::get<0>(src_image(row - radius, col - radius));
            std::get<1>(mirrored(row, col)) = std::get<1>(src_image(row - radius, col - radius));
            std::get<2>(mirrored(row, col)) = std::get<2>(src_image(row - radius, col - radius));
        }
    }
    for (uint row = 0; row < radius; ++row) {
        for (uint col = 0; col < src_image.n_cols; ++col) {
            std::get<0>(mirrored(radius - row - 1, col + radius)) = std::get<0>(mirrored(radius + 
                row, col + radius));
            std::get<1>(mirrored(radius - row - 1, col + radius)) = std::get<1>(mirrored(radius + 
                row, col + radius));
            std::get<2>(mirrored(radius - row - 1, col + radius)) = std::get<2>(mirrored(radius + 
                row, col + radius));
            std::get<0>(mirrored(mirrored.n_rows - radius + row, col + radius)) = std::get<0>(
                mirrored(mirrored.n_rows - radius - row - 1, col + radius));
            std::get<1>(mirrored(mirrored.n_rows - radius + row, col + radius)) = std::get<1>(
                mirrored(mirrored.n_rows - radius - row - 1, col + radius));
            std::get<2>(mirrored(mirrored.n_rows - radius + row, col + radius)) = std::get<2>(
                mirrored(mirrored.n_rows - radius - row - 1, col + radius));
        }
    }
    for (uint col = 0; col < radius; ++col) {
        for (uint row = 0; row < mirrored.n_rows; ++row) {
            std::get<0>(mirrored(row, radius - col - 1)) = std::get<0>(mirrored(row, 
                radius + col));
            std::get<1>(mirrored(row, radius - col - 1)) = std::get<1>(mirrored(row, 
                radius + col));
            std::get<2>(mirrored(row, radius - col - 1)) = std::get<2>(mirrored(row, 
                radius + col));
            std::get<0>(mirrored(row, mirrored.n_cols - radius + col)) = std::get<0>(mirrored(
                row, mirrored.n_cols - radius - col - 1));
            std::get<1>(mirrored(row, mirrored.n_cols - radius + col)) = std::get<1>(mirrored(
                row, mirrored.n_cols - radius - col - 1));
            std::get<2>(mirrored(row, mirrored.n_cols - radius + col)) = std::get<2>(mirrored(
                row, mirrored.n_cols - radius - col - 1));
        }
    }
    return mirrored;
}

Image unmirroring(Image &src_image, uint radius) {
    Image unmirrored(src_image.n_rows - 2 * radius, src_image.n_cols - 2 * radius);
    for (uint row = 0; row < unmirrored.n_rows; ++row) {
        for (uint col = 0; col < unmirrored.n_cols; ++col) {
            std::get<0>(unmirrored(row, col)) = std::get<0>(src_image(row + radius, col + 
                radius));
            std::get<1>(unmirrored(row, col)) = std::get<1>(src_image(row + radius, col + 
                radius));
            std::get<2>(unmirrored(row, col)) = std::get<2>(src_image(row + radius, col + 
                radius));
        }
    }
    return src_image;
}

Matrix<std::tuple<double, double, double>> RGBtoYCrCb(Image src_image) {
    Matrix<std::tuple<double, double, double>> result(src_image.n_rows, src_image.n_cols);
    uint r = 0, g = 0, b = 0;
    for (uint row = 0; row < result.n_rows; ++row) {
        for (uint col = 0; col < result.n_cols; ++col) {
            std::tie(r, g, b) = src_image(row, col);
            double dr = static_cast<double>(r) / 255.0;
            double dg = static_cast<double>(g) / 255.0;
            double db = static_cast<double>(b) / 255.0;
            double Y = 0.2989 * dr + 0.5866 * dg + 0.1145 * db;
            double Cr = 0.5000 * dr - 0.4184 * dg - 0.0816 * db;
            double Cb = -0.1687 * dr - 0.3313 * dg + 0.5000 * db;
            result(row, col) = std::make_tuple(Y, Cr, Cb);
        }
    }
    return result;
}

Image histogram_equalization(Image src_image) {
    std::vector<std::vector<int>> histogram(3);
    for (uint i = 0; i < 3; ++i) {
        histogram[i].resize(256);
    }
    for (uint row = 0; row < src_image.n_rows; ++row) {
        for (uint col = 0; col < src_image.n_cols; ++col) {
            ++histogram[0][std::get<0>(src_image(row, col))];
            ++histogram[1][std::get<1>(src_image(row, col))];
            ++histogram[2][std::get<2>(src_image(row, col))];
        }
    }
    std::vector<uint> mins(3, 0), maxs(3, 255);
    for (uint i = 0; i < 255; ++i) {
        for (uint j = 0; j < 3; ++j) {
            if (histogram[j][mins[j]] == 0) {
                ++mins[j];
            }
            if (histogram[j][maxs[j]] == 0) {
                --maxs[j];
            }
        }
    }
    Image result = src_image.deep_copy();
    std::vector<double> contrasts(3);
    for (uint i = 0; i < 3; ++i) {
        contrasts[i] = 255.0 / static_cast<double>(maxs[i] - mins[i]);
    }
    for (uint row = 0; row < result.n_rows; ++row) {
        for (uint col = 0; col < result.n_cols; ++col) {
            int r = contrasts[0] * (static_cast<int>(std::get<0>(result(row, col))) - 
                static_cast<int>(mins[0]));
            int g = contrasts[1] * (static_cast<int>(std::get<1>(result(row, col))) - 
                static_cast<int>(mins[1]));
            int b = contrasts[2] * (static_cast<int>(std::get<2>(result(row, col))) - 
                static_cast<int>(mins[2]));
            std::get<0>(result(row, col)) = std::min(std::max(r, 0), 255);
            std::get<1>(result(row, col)) = std::min(std::max(g, 0), 255);
            std::get<2>(result(row, col)) = std::min(std::max(b, 0), 255);
        }
    }
    return result;
}

Image white_balance(Image src_image) {
    Image result = histogram_equalization(src_image);
    Matrix<std::tuple<double, double, double>> YCrCb = RGBtoYCrCb(result);// std::tuple<Y, Cr, Cb>
    std::vector<std::tuple<uint, uint, double, double, double>> probable_wp;
    double sumY = 0, sumCr = 0, sumCb = 0;
    for (uint row = 0; row < YCrCb.n_rows; ++row) {
        for (uint col = 0; col < YCrCb.n_cols; ++col) {
            if ((std::get<0>(YCrCb(row, col)) > 210.0) && (std::get<1>(YCrCb(row, col)) >= -3) && 
                (std::get<2>(YCrCb(row, col)) <= 3)) {
                probable_wp.push_back(std::make_tuple(row, col, std::get<0>(YCrCb(row, col)), 
                    std::get<1>(YCrCb(row, col)), std::get<2>(YCrCb(row, col))));
                sumY += std::get<0>(YCrCb(row, col));
                sumCr += std::get<1>(YCrCb(row, col));
                sumCb += std::get<2>(YCrCb(row, col));
            }
        }
    }
    if (probable_wp.size() == 0) {
        return src_image;
    }
    std::vector<double> YCrCb_avg = {sumY / probable_wp.size(), sumCr / probable_wp.size(), 
        sumCb / probable_wp.size()};
    std::sort(probable_wp.begin(), probable_wp.end(), [&](std::tuple<uint, uint, double, 
        double, double> &a, std::tuple<uint, uint, double, double, double> &b){return 
        std::get<2>(a) > std::get<2>(b);});
    double norm = std::get<3>(probable_wp[0]) * std::get<3>(probable_wp[0]) + 
    std::get<4>(probable_wp[0]) * std::get<4>(probable_wp[0]);
    uint bright_index = 0;
    for (uint i = 0; i < probable_wp.size(); ++i) {
        double cur_norm = std::get<3>(probable_wp[i]) * std::get<3>(probable_wp[i]) + 
        std::get<4>(probable_wp[i]) * std::get<4>(probable_wp[i]);
        if (cur_norm < norm) {
            norm = cur_norm;
            bright_index = i;
        }
    }
    std::vector<double> YCrCb_bright = {std::get<2>(probable_wp[bright_index]), 
        std::get<3>(probable_wp[bright_index]), std::get<4>(probable_wp[bright_index])};
    std::vector<double> limits = {std::min(YCrCb_bright[0], YCrCb_avg[0]), 
        std::max(YCrCb_bright[0], YCrCb_avg[0]), std::min(YCrCb_bright[1], YCrCb_avg[1]), 
        std::max(YCrCb_bright[1], YCrCb_avg[1]), std::min(YCrCb_bright[2], YCrCb_avg[2]), 
        std::max(YCrCb_bright[2], YCrCb_avg[2])};
    std::vector<std::tuple<uint, uint>> whiteRGB;
    unsigned long long sum_red = 0, sum_green = 0, sum_blue = 0;
    for (uint row = 0; row < YCrCb.n_rows; ++row) {
        for (uint col = 0; col < YCrCb.n_cols; ++col) {
            double Y = std::get<0>(YCrCb(row, col));
            double Cr = std::get<1>(YCrCb(row, col));
            double Cb = std::get<2>(YCrCb(row, col));
            if ((limits[0] <= Y) && (Y <= limits[1]) && (limits[2] <= Cr) && (Cr <= limits[3]) && 
                (limits[4] <= Cb) && (Cb <= limits[5])) {
                whiteRGB.push_back(std::make_tuple(row, col));
                sum_red += std::get<0>(result(row, col));
                sum_green += std::get<1>(result(row, col));
                sum_blue += std::get<2>(result(row, col));
            }
        }
    }
    if (whiteRGB.size() == 0) {
        return src_image;
    }
    std::vector<double> W = {static_cast<double>(sum_red) / whiteRGB.size(), 
        static_cast<double>(sum_green) / whiteRGB.size(), 
        static_cast<double>(sum_blue) / whiteRGB.size()};
    unsigned long long sr = 0, sg = 0, sb = 0, y = 0;
    for (uint row = 0; row < result.n_rows; ++row) {
        for (uint col = 0; col < result.n_cols; ++col) {
            sr += std::get<0>(result(row, col));
            sg += std::get<1>(result(row, col));
            sb += std::get<2>(result(row, col));
            y += 0.299 * std::get<0>(result(row, col)) + 0.587 * std::get<1>(result(row, col)) + 
            0.114 * std::get<2>(result(row, col));
        }
    }
    sr /= static_cast<double>(result.n_rows * result.n_cols);
    sg /= static_cast<double>(result.n_rows * result.n_cols);
    sb /= static_cast<double>(result.n_rows * result.n_cols);
    y /= static_cast<double>(result.n_rows * result.n_cols);
    double Yw = 0.299 * W[0] + 0.587 * W[1] + 0.114 * W[2];
    std::vector<double> ScaleFactors = {static_cast<double>(Yw) / W[0], 
        static_cast<double>(Yw) / W[1], static_cast<double>(Yw) / W[2], 
        static_cast<double>(y) / sr, static_cast<double>(y) / sg, static_cast<double>(y) / sb};
    uint R = std::min(std::max(YCrCb_avg[0] + 1.4022 * YCrCb_avg[1], 0.0), 255.0);
    uint G = std::min(std::max(YCrCb_avg[0] - 0.7145 * YCrCb_avg[1] - 0.3456 * YCrCb_avg[2], 0.0), 
        255.0);
    uint B = std::min(std::max(YCrCb_avg[0] + 1.7710 * YCrCb_avg[2], 0.0), 255.0);
    double Rfactor = ScaleFactors[0];
    double Gfactor = ScaleFactors[1];
    double Bfactor = ScaleFactors[2];
    int color_cast = 0;
    if ((B + 3 >= G) && (B > R)) {
        color_cast = 3;
        Bfactor = ScaleFactors[5];
    } else {
        if ((G + 3 > R) && (R > B)) {
            color_cast = 2;
            Gfactor = ScaleFactors[4];
        } else {
            if ((R > G) && (G > B)) {
                color_cast = 1;
                Rfactor = ScaleFactors[3];
            }
        }
    }
    if (color_cast > 0) {
        for (uint row = 0; row < result.n_rows; ++row) {
            for (uint col = 0; col < result.n_cols; ++col) {
                uint r = std::get<0>(result(row, col));
                uint g = std::get<1>(result(row, col));
                uint b = std::get<2>(result(row, col));
                r *= Rfactor;
                g *= Gfactor;
                b *= Bfactor;
                std::get<0>(result(row, col)) = std::min(r, 255u);
                std::get<1>(result(row, col)) = std::min(g, 255u);
                std::get<2>(result(row, col)) = std::min(b, 255u);
            }
        }
    }
    return result;
}

Image align(Image srcImage, bool isPostprocessing, std::string postprocessingType, 
    double fraction, bool isMirror, 
            bool isInterp, bool isSubpixel, double subScale)
{
    int threshold1 = 10, threshold2 = 30;
    Image edge_map = canny(srcImage, threshold1, threshold2, isMirror);
    double cut_scale = 0.1;
    RGB channels = cut_image(srcImage, edge_map, cut_scale);
    int delta_shift = 0.05 * srcImage.n_cols;
    Image RG = find_mse(channels.R, channels.G, delta_shift, 2);
    Image BGR = find_mse(channels.B, RG, delta_shift, 0);
    if (isPostprocessing) {
        if (postprocessingType == "--gray-world") {
            BGR = gray_world(BGR);
        }
        if (postprocessingType == "--unsharp") {
            BGR = unsharp(BGR);
        }
        if (postprocessingType == "--autocontrast") {
            BGR = autocontrast(BGR, fraction);
        }
        if (postprocessingType == "--white-balance") {
            BGR = white_balance(BGR);
        }
    }
    return BGR;
}

Image sobel_x(Image src_image) {
    Matrix<double> kernel = {{-1, 0, 1},
                             {-2, 0, 2},
                             {-1, 0, 1}};
    return custom(src_image, kernel);
}

Image sobel_y(Image src_image) {
    Matrix<double> kernel = {{ 1,  2,  1},
                             { 0,  0,  0},
                             {-1, -2, -1}};
    return custom(src_image, kernel);
}

Image unsharp(Image src_image) {
    Matrix<double> kernel = {{-1.0 / 6.0, -2.0 / 3.0, -1.0 / 6.0}, 
                             {-2.0 / 3.0,  13.0 / 3.0, -2.0 / 3.0}, 
                             {-1.0 / 6.0, -2.0 / 3.0, -1.0 / 6.0}};
    //Matrix<double> kernel = {{-0.1, -0.2, -0.1}, 
    //                         {-0.2,  2.2, -0.2}, 
    //                         {-0.1, -0.2, -0.1}};
    return custom(src_image, kernel);
}

Image gray_world(Image src_image) {
    Image result = src_image.deep_copy();
    unsigned long long sr = 0, sg = 0, sb = 0;
    double srgb = 0;
    for (uint row = 0; row < result.n_rows; ++row) {
        for (uint col = 0; col < result.n_cols; ++col) {
            sr += std::get<0>(result(row, col));
            sg += std::get<1>(result(row, col));
            sb += std::get<2>(result(row, col));
        }
    }
    sr /= static_cast<double>(result.n_rows * result.n_cols);
    sg /= static_cast<double>(result.n_rows * result.n_cols);
    sb /= static_cast<double>(result.n_rows * result.n_cols);
    srgb = (sr + sg + sb) / 3.0;
    uint r = 0, g = 0, b = 0;
    for (uint row = 0; row < result.n_rows; ++row) {
        for (uint col = 0; col < result.n_cols; ++col) {
            r = std::get<0>(result(row, col));
            g = std::get<1>(result(row, col));
            b = std::get<2>(result(row, col));
            r *= srgb / sr;
            g *= srgb / sg;
            b *= srgb / sb;
            std::get<0>(result(row, col)) = std::min(r, 255u);
            std::get<1>(result(row, col)) = std::min(g, 255u);
            std::get<2>(result(row, col)) = std::min(b, 255u);
        }
    }
    return result;
}

Image resize(Image src_image, double scale) {
    return src_image;
}

template<typename ReturnT, typename KernelT>//based on class BoxxFilterOp
class KernelOp
{
    Matrix<KernelT> kernel;
public:
    uint radius;
    KernelOp(Matrix<KernelT> k, uint r = 1): kernel(k), radius(r) {}
    tuple<ReturnT, ReturnT, ReturnT> operator () (const Image &m) const
    {
        uint size = 2 * radius + 1;
        ReturnT r = 0, g = 0, b = 0;
        KernelT sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;
        for (uint i = 0; i < size; ++i) {
            for (uint j = 0; j < size; ++j) {
                // Tie is useful for taking elements from tuple
                tie(r, g, b) = m(i, j);
                sum_r += r * kernel(i, j);
                sum_g += g * kernel(i, j);
                sum_b += b * kernel(i, j);
            }
        }
        if (typeid(sum_r) != typeid(r)) {
            r = std::min(std::max(static_cast<int>(sum_r), 0), 255);
            g = std::min(std::max(static_cast<int>(sum_g), 0), 255);
            b = std::min(std::max(static_cast<int>(sum_b), 0), 255);
        } else {
            r = sum_r;
            g = sum_g;
            b = sum_b;
        }
        return make_tuple(r, g, b);
    }
};

Image custom(Image src_image, Matrix<double> kernel) {
    // Function custom is useful for making concrete linear filtrations
    // like gaussian or sobel. So, we assume that you implement custom
    // and then implement other filtrations using this function.
    // sobel_x and sobel_y are given as an example.
    Image mirrored = mirroring(src_image, kernel.n_rows / 2);
    Image result = mirrored.unary_map(KernelOp<uint, double>(kernel, kernel.n_rows / 2));
    //kernel should have square shape
    return unmirroring(result, kernel.n_rows / 2);
}

Image autocontrast(Image src_image, double fraction) {
    std::vector<uint> histogram(256, 0);
    for (uint row = 0; row < src_image.n_rows; ++row) {
        for (uint col = 0; col < src_image.n_cols; ++col) {
            uint r = 0.2125 * std::get<0>(src_image(row, col));
            uint g = 0.7154 * std::get<1>(src_image(row, col));
            uint b = 0.0721 * std::get<2>(src_image(row, col));
            ++histogram[std::min(r + g + b, 255u)];
        }
    }
    uint low_sum = src_image.n_rows * src_image.n_cols * fraction;
    uint cumsum = 0, ymin = 0, ymax = 255;
    for (uint i = 0; i < 255; ++i) {
        if (histogram[ymin] == 0) {
            ++ymin;
        }
        if (histogram[ymax] == 0) {
            --ymax;
        }
    }
    if (fraction > 0) {
        while ((cumsum < low_sum) && (ymin < 255)) {
            cumsum += histogram[ymin++];
        }
        ymax = ymin;
        --ymin;
        low_sum = src_image.n_rows * src_image.n_cols - low_sum;
        while((cumsum < low_sum) && (ymax < 256)) {
            cumsum += histogram[ymax++];
        }
        --ymax;
    }
    Image result = src_image.deep_copy();
    double contrast = 255.0 / static_cast<double>(ymax - ymin);
    for (uint row = 0; row < result.n_rows; ++row) {
        for (uint col = 0; col < result.n_cols; ++col) {
            int r = contrast * (static_cast<int>(std::get<0>(result(row, col))) - 
                static_cast<int>(ymin));
            int g = contrast * (static_cast<int>(std::get<1>(result(row, col))) - 
                static_cast<int>(ymin));
            int b = contrast * (static_cast<int>(std::get<2>(result(row, col))) - 
                static_cast<int>(ymin));
            std::get<0>(result(row, col)) = std::min(std::max(r, 0), 255);
            std::get<1>(result(row, col)) = std::min(std::max(g, 0), 255);
            std::get<2>(result(row, col)) = std::min(std::max(b, 0), 255);
        }
    }
    return result;
}

double gauss_gen_2d(double x, double y, double sigma, double pi) {
    return 1.0 / (2.0 * pi * sigma * sigma) * std::exp(-(x * x + y * y) / (2 * sigma * sigma));
}

Image gaussian(Image src_image, double sigma, int radius)  {
    double pi = 2.0 * std::asin(1.0);
    int size = 2 * radius + 1;
    Matrix<double> kernel(size, size);
    double sum = 0.0;
    for (int i = -radius; i < radius + 1; ++i) {
        for (int j = -radius; j < radius + 1; ++j) {
            double num = gauss_gen_2d(i, j, sigma, pi);
            kernel(radius + i, radius + j) = num;
            sum += num;
        }
    }
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            kernel(i, j) /= sum;
        }
    }
    return custom(src_image, kernel);
}

double gauss_gen_1d(double x, double sigma, double pi) {
    return 1.0 / (std::sqrt(2.0 * pi) * sigma) * std::exp(-(x * x) / (2 * sigma * sigma));
}

Image gaussian_separable(Image src_image, double sigma, int radius, bool mirrorfunc) {
    double pi = 2.0 * std::asin(1.0);
    int size = 2 * radius + 1;
    Matrix<double> kernel(1, size);
    double sum = 0.0;
    for (int i = -radius; i < radius + 1; ++i) {
        double num = gauss_gen_1d(i, sigma, pi);
        kernel(0, radius + i) = num;
        sum += num;
    }
    for (int i = 0; i < size; ++i) {
        kernel(0, i) /= sum;
    }
    Image mirrored = src_image;
    if (mirrorfunc) {
        mirrored = mirroring(src_image, radius);
    }
    Image tmp = mirrored.deep_copy();
    uint end_col = tmp.n_cols - radius;
    for (uint row = 0; row < tmp.n_rows; ++row) {
        for (uint col = radius; col < end_col; ++col) {
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;
            uint r = 0, g = 0, b = 0;
            for (int i = -radius; i < radius + 1; ++i) {
                tie(r, g, b) = mirrored(row, col + i);
                sum_r += r * kernel(0, radius + i);
                sum_g += g * kernel(0, radius + i);
                sum_b += b * kernel(0, radius + i);
            }
            r = sum_r;
            g = sum_g;
            b = sum_b;
            tmp(row, col) = make_tuple(r, g, b);
        }
    }
    Image result = tmp.deep_copy();
    uint end_row = result.n_rows - radius;
    for (uint row = radius; row < end_row; ++row) {
        for (uint col = 0; col < result.n_cols; ++col) {
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;
            uint r = 0, g = 0, b = 0;
            for (int i = -radius; i < radius + 1; ++i) {
                tie(r, g, b) = tmp(row + i, col);
                sum_r += r * kernel(0, radius + i);
                sum_g += g * kernel(0, radius + i);
                sum_b += b * kernel(0, radius + i);
            }
            r = sum_r;
            g = sum_g;
            b = sum_b;
            result(row, col) = make_tuple(r, g, b);
        }
    }
    if (mirrorfunc) {
        return unmirroring(result, radius);
    } else {
        return result;
    }
}

Image median(Image src_image, int radius) {
    Image result = mirroring(src_image, radius);
    Image tmp = result.deep_copy();
    std::vector<std::vector<uint>> kernel(3);
    uint size = 2 * radius + 1;
    uint placed = size * size;
    for (uint i = 0; i < 3; ++i) {
        kernel[i].resize(placed);
    }
    uint end_row = result.n_rows - radius;
    uint end_col = result.n_cols - radius;
    for (uint row = radius; row < end_row; ++row) {
        for (uint col = radius; col < end_col; ++col) {
            uint i = row - radius;
            uint j = col - radius;
            for (uint elem = 0; elem < placed; ++elem) {
                kernel[0][elem] = std::get<0>(tmp(i, j));
                kernel[1][elem] = std::get<1>(tmp(i, j));
                kernel[2][elem] = std::get<2>(tmp(i, j));
                ++j;
                if ((elem + 1) % size == 0) {
                    ++i;
                    j -= size; 
                }
            }
            std::sort(kernel[0].begin(), kernel[0].end());
            std::sort(kernel[1].begin(), kernel[1].end());
            std::sort(kernel[2].begin(), kernel[2].end());
            std::get<0>(result(row, col)) = kernel[0][kernel[0].size() / 2];
            std::get<1>(result(row, col)) = kernel[1][kernel[1].size() / 2];
            std::get<2>(result(row, col)) = kernel[2][kernel[2].size() / 2];
        }
    }
    return unmirroring(result, radius);
}

Image median_linear(Image src_image, int radius) {
    Image result = mirroring(src_image, radius);
    Image tmp = result.deep_copy();
    int end_row = result.n_rows - radius - 1;
    int end_col = result.n_cols - radius;
    std::vector<std::vector<int>> histogram(3);
    for (int i = 0; i < 3; ++i) {
        histogram[i].resize(256);
    }
    int row = radius, col = radius, size = 2 * radius + 1;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            ++histogram[0][std::get<0>(tmp(i, j))];
            ++histogram[1][std::get<1>(tmp(i, j))];
            ++histogram[2][std::get<2>(tmp(i, j))];
        }
    }
    int stop_num = size * size / 2;
    std::vector<int> cumsum(3, 0), median(3, 0);
    while (cumsum[0] < stop_num) {
        cumsum[0] += histogram[0][median[0]++];
    }
    --median[0];
    while (cumsum[1] < stop_num) {
        cumsum[1] += histogram[1][median[1]++];
    }
    --median[1];
    while (cumsum[2] < stop_num) {
        cumsum[2] += histogram[2][median[2]++];
    }
    --median[2];
    std::get<0>(result(row, col)) = median[0];
    std::get<1>(result(row, col)) = median[1];
    std::get<2>(result(row, col)) = median[2];
    bool isright = true, analysed = true, can_move = true;
    while (true) {
        cumsum[0] = 0;
        cumsum[1] = 0;
        cumsum[2] = 0;
        median[0] = 0;
        median[1] = 0;
        median[2] = 0;
        analysed = false;
        if (isright && can_move) {
            for (int i = row - radius; i < row + radius + 1; ++i) {
                --histogram[0][std::get<0>(tmp(i, col - radius))];
                --histogram[1][std::get<1>(tmp(i, col - radius))];
                --histogram[2][std::get<2>(tmp(i, col - radius))];
            }
            ++col;
            for (int i = row - radius; i < row + radius + 1; ++i) {
                ++histogram[0][std::get<0>(tmp(i, col + radius))];
                ++histogram[1][std::get<1>(tmp(i, col + radius))];
                ++histogram[2][std::get<2>(tmp(i, col + radius))];
            }
            analysed = true;
            if (col == end_col - 1) {
                can_move = false;
            }
        }
        if (!isright && can_move) {
            for (int i = row - radius; i < row + radius + 1; ++i) {
                --histogram[0][std::get<0>(tmp(i, col + radius))];
                --histogram[1][std::get<1>(tmp(i, col + radius))];
                --histogram[2][std::get<2>(tmp(i, col + radius))];
            }
            --col;
            for (int i = row - radius; i < row + radius + 1; ++i) {
                ++histogram[0][std::get<0>(tmp(i, col - radius))];
                ++histogram[1][std::get<1>(tmp(i, col - radius))];
                ++histogram[2][std::get<2>(tmp(i, col - radius))];
            }
            analysed = true;
            if (col == radius) {
                can_move = false;
            }
        }
        if (!analysed) {
            isright ^= true;
            for (int j = col - radius; j < col + radius + 1; ++j) {
                --histogram[0][std::get<0>(tmp(row - radius, j))];
                --histogram[1][std::get<1>(tmp(row - radius, j))];
                --histogram[2][std::get<2>(tmp(row - radius, j))];
            }
            if (row < end_row) {
                ++row;
            } else {
                break;
            }
            for (int j = col - radius; j < col + radius + 1; ++j) {
                ++histogram[0][std::get<0>(tmp(row + radius, j))];
                ++histogram[1][std::get<1>(tmp(row + radius, j))];
                ++histogram[2][std::get<2>(tmp(row + radius, j))];
            }
            can_move = true;
        }
        while (cumsum[0] < stop_num) {
            cumsum[0] += histogram[0][median[0]++];
        }
        --median[0];
        while (cumsum[1] < stop_num) {
            cumsum[1] += histogram[1][median[1]++];
        }
        --median[1];
        while (cumsum[2] < stop_num) {
            cumsum[2] += histogram[2][median[2]++];
        }
        --median[2];
        std::get<0>(result(row, col)) = median[0];
        std::get<1>(result(row, col)) = median[1];
        std::get<2>(result(row, col)) = median[2];
    }
    return unmirroring(result, radius);
}

Image median_const(Image src_image, int radius) {
    return src_image;
}

void nonmax_elimination(Matrix<tuple<double, double, bool>> &grad, double pi) {
    for (uint row = 0; row < grad.n_rows; ++row) {
        std::get<0>(grad(row, 0)) = 0;
        std::get<1>(grad(row, 0)) = 4 * pi;
        std::get<0>(grad(row, grad.n_cols - 1)) = 0;
        std::get<1>(grad(row, grad.n_cols - 1)) = 4 * pi;
    }
    for (uint col = 0; col < grad.n_cols; ++col) {
        std::get<0>(grad(0, col)) = 0;
        std::get<1>(grad(0, col)) = 4 * pi;
        std::get<0>(grad(grad.n_rows - 1, col)) = 0;
        std::get<1>(grad(grad.n_rows - 1, col)) = 4 * pi;
    }
    for (uint row = 1; row < grad.n_rows - 1; ++row) {
        for (uint col = 1; col < grad.n_cols - 1; ++col) {
            int a_row = 0, a_col = 0, b_row = 0, b_col = 0;
            double ahead = 0.0, behind = 0.0, norm = std::get<0>(grad(row, col));
            double angle = std::get<1>(grad(row, col));
            if (abs(angle) > pi) {
                continue;
            }
            if ((angle < (-7 * pi / 8)) || (angle > (7 * pi / 8))) {
                //left
                a_row = row;
                a_col = col - 1;
                b_row = row;
                b_col = col + 1;
                ahead = std::get<0>(grad(a_row, a_col));
                behind = std::get<0>(grad(b_row, b_col));
                if ((norm > ahead) && (norm > behind)) {
                    continue;
                } else {
                    std::get<0>(grad(row, col)) = 0;
                    std::get<1>(grad(row, col)) = 4 * pi;
                }
                continue;
            }
            if ((angle >= (-7 * pi / 8)) && (angle < (-5 * pi / 8))) {
                //left and down
                a_row = row + 1;
                a_col = col - 1;
                b_row = row - 1;
                b_col = col + 1;
                ahead = std::get<0>(grad(a_row, a_col));
                behind = std::get<0>(grad(b_row, b_col));
                if ((norm > ahead) && (norm > behind)) {
                    continue;
                } else {
                    std::get<0>(grad(row, col)) = 0;
                    std::get<1>(grad(row, col)) = 4 * pi;
                }
                continue;
            }
            if ((angle >= (-5 * pi / 8)) && (angle < (-3 * pi / 8))) {
                //down
                a_row = row + 1;
                a_col = col;
                b_row = row - 1;
                b_col = col;
                ahead = std::get<0>(grad(a_row, a_col));
                behind = std::get<0>(grad(b_row, b_col));
                if ((norm > ahead) && (norm > behind)) {
                    continue;
                } else {
                    std::get<0>(grad(row, col)) = 0;
                    std::get<1>(grad(row, col)) = 4 * pi;
                }
                continue;
            }
            if ((angle >= (-3 * pi / 8)) && (angle < (-1 * pi / 8))) {
                //down and right
                a_row = row + 1;
                a_col = col + 1;
                b_row = row - 1;
                b_col = col - 1;
                ahead = std::get<0>(grad(a_row, a_col));
                behind = std::get<0>(grad(b_row, b_col));
                if ((norm > ahead) && (norm > behind)) {
                    continue;
                } else {
                    std::get<0>(grad(row, col)) = 0;
                    std::get<1>(grad(row, col)) = 4 * pi;
                }
                continue;
            }
            if ((angle >= (-1 * pi / 8)) && (angle < (1 * pi / 8))) {
                //right
                a_row = row;
                a_col = col + 1;
                b_row = row;
                b_col = col - 1;
                ahead = std::get<0>(grad(a_row, a_col));
                behind = std::get<0>(grad(b_row, b_col));
                if ((norm > ahead) && (norm > behind)) {
                    continue;
                } else {
                    std::get<0>(grad(row, col)) = 0;
                    std::get<1>(grad(row, col)) = 4 * pi;
                }
                continue;
            }
            if ((angle >= (1 * pi / 8)) && (angle < (3 * pi / 8))) {
                //right and up
                a_row = row - 1;
                a_col = col + 1;
                b_row = row + 1;
                b_col = col - 1;
                ahead = std::get<0>(grad(a_row, a_col));
                behind = std::get<0>(grad(b_row, b_col));
                if ((norm > ahead) && (norm > behind)) {
                    continue;
                } else {
                    std::get<0>(grad(row, col)) = 0;
                    std::get<1>(grad(row, col)) = 4 * pi;
                }
                continue;
            }
            if ((angle >= (3 * pi / 8)) && (angle < (5 * pi / 8))) {
                //up
                a_row = row - 1;
                a_col = col;
                b_row = row + 1;
                b_col = col;
                ahead = std::get<0>(grad(a_row, a_col));
                behind = std::get<0>(grad(b_row, b_col));
                if ((norm > ahead) && (norm > behind)) {
                    continue;
                } else {
                    std::get<0>(grad(row, col)) = 0;
                    std::get<1>(grad(row, col)) = 4 * pi;
                }
                continue;
            }
            if ((angle >= (5 * pi / 8)) && (angle < (7 * pi / 8))) {
                //up and left
                a_row = row - 1;
                a_col = col - 1;
                b_row = row + 1;
                b_col = col + 1;
                ahead = std::get<0>(grad(a_row, a_col));
                behind = std::get<0>(grad(b_row, b_col));
                if ((norm > ahead) && (norm > behind)) {
                    continue;
                } else {
                    std::get<0>(grad(row, col)) = 0;
                    std::get<1>(grad(row, col)) = 4 * pi;
                }
                continue;
            }
        }
    }
}

struct Coord
{
    uint row, col;
    Coord(int r = 0, int c = 0): row(r), col(c) {}
};

std::vector<Coord> double_thresholding(Matrix<tuple<double, double, bool>> &grad, int threshold1, 
    int threshold2, double pi) {
    std::vector<Coord> strongs;
    for (uint row = 0; row < grad.n_rows; ++row) {
        for (uint col = 0; col < grad.n_cols; ++col) {
            if (static_cast<int>(std::get<0>(grad(row, col))) < threshold1) {
                std::get<0>(grad(row, col)) = 0;
                std::get<1>(grad(row, col)) = 4 * pi;
                continue;
            }
            if (static_cast<int>(std::get<0>(grad(row, col))) > threshold2) {
                std::get<2>(grad(row, col)) = true;
                strongs.push_back(Coord(row, col));
                continue;
            }
        }
    }
    return strongs;
}

void dfs(Matrix<tuple<double, double, bool>> &grad, const Coord &position, double pi) {
    int a_row = 0, a_col = 0, b_row = 0, b_col = 0;
    double angle = std::get<1>(grad(position.row, position.col)) + pi / 2.0;
    if (angle > pi) {
        angle -= 2 * pi;
    }
    if (abs(angle) > pi) {
        return;
    }
    if ((angle < (-7 * pi / 8)) || (angle > (7 * pi / 8))) {
        //left
        a_row = position.row;
        a_col = position.col - 1;
        b_row = position.row;
        b_col = position.col + 1;
        if ((a_col < 0) || (b_col >= static_cast<int>(grad.n_cols))) {
            return;
        }
        if ((std::get<0>(grad(a_row, a_col)) > 0) && (std::get<2>(grad(a_row, a_col)) == false)) {
            std::get<2>(grad(a_row, a_col)) = true;
            dfs(grad, Coord(a_row, a_col), pi);
        }
        if ((std::get<0>(grad(b_row, b_col)) > 0) && (std::get<2>(grad(b_row, b_col)) == false)) {
            std::get<2>(grad(b_row, b_col)) = true;
            dfs(grad, Coord(b_row, b_col), pi);
        }
        return;
    }
    if ((angle >= (-7 * pi / 8)) && (angle < (-5 * pi / 8))) {
        //left and down
        a_row = position.row + 1;
        a_col = position.col - 1;
        b_row = position.row - 1;
        b_col = position.col + 1;
        if ((a_row >= static_cast<int>(grad.n_rows)) || (a_col < 0) || (b_row < 0) || (b_col >= 
            static_cast<int>(grad.n_cols))) {
            return;
        }
        if ((std::get<0>(grad(a_row, a_col)) > 0) && (std::get<2>(grad(a_row, a_col)) == false)) {
            std::get<2>(grad(a_row, a_col)) = true;
            dfs(grad, Coord(a_row, a_col), pi);
        }
        if ((std::get<0>(grad(b_row, b_col)) > 0) && (std::get<2>(grad(b_row, b_col)) == false)) {
            std::get<2>(grad(b_row, b_col)) = true;
            dfs(grad, Coord(b_row, b_col), pi);
        }
        return;
    }
    if ((angle >= (-5 * pi / 8)) && (angle < (-3 * pi / 8))) {
        //down
        a_row = position.row + 1;
        a_col = position.col;
        b_row = position.row - 1;
        b_col = position.col;
        if ((a_row >= static_cast<int>(grad.n_rows)) || (b_row < 0)) {
            return;
        }
        if ((std::get<0>(grad(a_row, a_col)) > 0) && (std::get<2>(grad(a_row, a_col)) == false)) {
            std::get<2>(grad(a_row, a_col)) = true;
            dfs(grad, Coord(a_row, a_col), pi);
        }
        if ((std::get<0>(grad(b_row, b_col)) > 0) && (std::get<2>(grad(b_row, b_col)) == false)) {
            std::get<2>(grad(b_row, b_col)) = true;
            dfs(grad, Coord(b_row, b_col), pi);
        }
        return;
    }
    if ((angle >= (-3 * pi / 8)) && (angle < (-1 * pi / 8))) {
        //down and right
        a_row = position.row + 1;
        a_col = position.col + 1;
        b_row = position.row - 1;
        b_col = position.col - 1;
        if ((a_row >= static_cast<int>(grad.n_rows)) || (a_col >= 
            static_cast<int>(grad.n_cols)) || (b_row < 0) || (b_col < 0)) {
            return;
        }
        if ((std::get<0>(grad(a_row, a_col)) > 0) && (std::get<2>(grad(a_row, a_col)) == false)) {
            std::get<2>(grad(a_row, a_col)) = true;
            dfs(grad, Coord(a_row, a_col), pi);
        }
        if ((std::get<0>(grad(b_row, b_col)) > 0) && (std::get<2>(grad(b_row, b_col)) == false)) {
            std::get<2>(grad(b_row, b_col)) = true;
            dfs(grad, Coord(b_row, b_col), pi);
        }
        return;
    }
    if ((angle >= (-1 * pi / 8)) && (angle < (1 * pi / 8))) {
        //right
        a_row = position.row;
        a_col = position.col + 1;
        b_row = position.row;
        b_col = position.col - 1;
        if ((a_col >= static_cast<int>(grad.n_cols)) || (b_col < 0)) {
            return;
        }
        if ((std::get<0>(grad(a_row, a_col)) > 0) && (std::get<2>(grad(a_row, a_col)) == false)) {
            std::get<2>(grad(a_row, a_col)) = true;
            dfs(grad, Coord(a_row, a_col), pi);
        }
        if ((std::get<0>(grad(b_row, b_col)) > 0) && (std::get<2>(grad(b_row, b_col)) == false)) {
            std::get<2>(grad(b_row, b_col)) = true;
            dfs(grad, Coord(b_row, b_col), pi);
        }
        return;
    }
    if ((angle >= (1 * pi / 8)) && (angle < (3 * pi / 8))) {
        //right and up
        a_row = position.row - 1;
        a_col = position.col + 1;
        b_row = position.row + 1;
        b_col = position.col - 1;
        if ((a_row < 0) || (a_col >= static_cast<int>(grad.n_cols)) || 
            (b_row >= static_cast<int>(grad.n_rows)) || (b_col < 0)) {
            return;
        }
        if ((std::get<0>(grad(a_row, a_col)) > 0) && (std::get<2>(grad(a_row, a_col)) == false)) {
            std::get<2>(grad(a_row, a_col)) = true;
            dfs(grad, Coord(a_row, a_col), pi);
        }
        if ((std::get<0>(grad(b_row, b_col)) > 0) && (std::get<2>(grad(b_row, b_col)) == false)) {
            std::get<2>(grad(b_row, b_col)) = true;
            dfs(grad, Coord(b_row, b_col), pi);
        }
        return;
    }
    if ((angle >= (3 * pi / 8)) && (angle < (5 * pi / 8))) {
        //up
        a_row = position.row - 1;
        a_col = position.col;
        b_row = position.row + 1;
        b_col = position.col;
        if ((a_row < 0) || (b_row >= static_cast<int>(grad.n_rows))) {
            return;
        }
        if ((std::get<0>(grad(a_row, a_col)) > 0) && (std::get<2>(grad(a_row, a_col)) == false)) {
            std::get<2>(grad(a_row, a_col)) = true;
            dfs(grad, Coord(a_row, a_col), pi);
        }
        if ((std::get<0>(grad(b_row, b_col)) > 0) && (std::get<2>(grad(b_row, b_col)) == false)) {
            std::get<2>(grad(b_row, b_col)) = true;
            dfs(grad, Coord(b_row, b_col), pi);
        }
        return;
    }
    if ((angle >= (5 * pi / 8)) && (angle < (7 * pi / 8))) {
        //up and left
        a_row = position.row - 1;
        a_col = position.col - 1;
        b_row = position.row + 1;
        b_col = position.col + 1;
        if ((a_row < 0) || (a_col < 0) || (b_row >= static_cast<int>(grad.n_rows)) || 
            (b_col >= static_cast<int>(grad.n_cols))) {
            return;
        }
        if ((std::get<0>(grad(a_row, a_col)) > 0) && (std::get<2>(grad(a_row, a_col)) == false)) {
            std::get<2>(grad(a_row, a_col)) = true;
            dfs(grad, Coord(a_row, a_col), pi);
        }
        if ((std::get<0>(grad(b_row, b_col)) > 0) && (std::get<2>(grad(b_row, b_col)) == false)) {
            std::get<2>(grad(b_row, b_col)) = true;
            dfs(grad, Coord(b_row, b_col), pi);
        }
        return;
    }
}

void histeresis(Matrix<tuple<double, double, bool>> &grad, std::vector<Coord> &strongs, 
    double pi) {
    for (auto strong = strongs.cbegin(); strong != strongs.cend(); ++strong) {
        dfs(grad, *strong, pi);
    }
    for (uint row = 0; row < grad.n_rows; ++row) {
        for (uint col = 0; col < grad.n_cols; ++col) {
            if (std::get<2>(grad(row, col)) == false) {
                std::get<0>(grad(row, col)) = 0;
                std::get<1>(grad(row, col)) = 4 * pi;
            }
        }
    }
}

Image canny(Image src_image, int threshold1, int threshold2, bool mirrorfunc) {
    Image smoothered = gaussian_separable(src_image, 1.4, 2, mirrorfunc);
    if (mirrorfunc) {
        smoothered = mirroring(smoothered, 1);
    }
    Matrix<tuple<int, int, int>> sobelX = smoothered.unary_map(KernelOp<int, int>({{-1, 0, 1}, 
        {-2, 0, 2}, {-1, 0, 1}}, 1));
    Matrix<tuple<int, int, int>> sobelY = smoothered.unary_map(KernelOp<int, int>({{1, 2, 1}, 
        {0, 0, 0}, {-1, -2, -1}}, 1));
    if (mirrorfunc) {
        sobelX = sobelX.submatrix(1, 1, sobelX.n_rows - 2, sobelX.n_cols - 2);
        sobelY = sobelY.submatrix(1, 1, sobelY.n_rows - 2, sobelY.n_cols - 2);
    }
    Matrix<tuple<double, double, bool>> gradient(sobelX.n_rows, sobelX.n_cols);
    for (uint row = 0; row < gradient.n_rows; ++row) {
        for (uint col = 0; col < gradient.n_cols; ++col) {
            int sobel_x_pixel = std::get<0>(sobelX(row, col)), sobel_y_pixel = 
            std::get<0>(sobelY(row, col));
            std::get<0>(gradient(row, col)) = std::sqrt(sobel_x_pixel * sobel_x_pixel + 
                sobel_y_pixel * sobel_y_pixel);
            std::get<1>(gradient(row, col)) = std::atan2(sobel_y_pixel, sobel_x_pixel);
            std::get<2>(gradient(row, col)) = false;
        }
    }
    double pi = 2.0 * std::asin(1.0);
    nonmax_elimination(gradient, pi);
    std::vector<Coord> strongs = double_thresholding(gradient, threshold1, threshold2, pi);
    histeresis(gradient, strongs, pi);
    Image result(gradient.n_rows, gradient.n_cols);
    for (uint row = 0; row < result.n_rows; ++row) {
        for (uint col = 0; col < result.n_cols; ++col) {
            uint tmp = std::abs(std::get<0>(gradient(row, col)));
            tmp = std::min(tmp, 255u);
            std::get<0>(result(row, col)) = tmp;
            std::get<1>(result(row, col)) = tmp;
            std::get<2>(result(row, col)) = tmp;
        }
    }
    return result;
}
