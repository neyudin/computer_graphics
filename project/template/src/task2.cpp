#include <string>
#include <vector>
#include <fstream>
#include <cassert>
#include <iostream>
#include <cmath>

#include "classifier.h"
#include "EasyBMP.h"
#include "linear.h"
#include "argvparser.h"
#include "matrix.h"

using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;

using CommandLineProcessing::ArgvParser;

float PI = 2 * std::asin(1.0);

typedef vector<pair<BMP*, int> > TDataSet;
typedef vector<pair<string, int> > TFileList;
typedef vector<pair<vector<float>, int> > TFeatures;

Matrix<std::tuple<uint, uint, uint>> fromBMPtoMatrix(BMP *src_image) {
    Matrix<std::tuple<uint, uint, uint>> result(src_image->TellHeight(), src_image->TellWidth());
    for (uint i = 0; i < result.n_rows; ++i) {
        for (uint j = 0; j < result.n_cols; ++j) {
            RGBApixel p = src_image->GetPixel(j, i);
            result(i, j) = std::make_tuple(p.Red, p.Green, p.Blue);
        }
    }
    return result;
}

Matrix<std::tuple<uint, uint, uint>> grayscale(Matrix<std::tuple<uint, uint, uint>> &src_image) {
    Matrix<std::tuple<uint, uint, uint>> result(src_image.n_rows, src_image.n_cols);
    for (uint row = 0; row < result.n_rows; ++row) {
        for (uint col = 0; col < result.n_cols; ++col) {
            uint r = 0, g = 0, b = 0;
            std::tie(r, g, b) = src_image(row, col);
            uint intensity = 0.299 * r + 0.587 * g + 0.114 * b;
            result(row, col) = std::make_tuple(intensity, intensity, intensity);
        }
    }
    return result;
}

template<typename ReturnT, typename PixelT, typename KernelT>//based on class BoxxFilterOp
class KernelOp
{
    Matrix<KernelT> kernel;
public:
    uint vert_radius, hor_radius;
    KernelOp(Matrix<KernelT> k, uint vr = 1, uint hr = 1): kernel(k), vert_radius(vr),
        hor_radius(hr) {}
    std::tuple<ReturnT, ReturnT, ReturnT> operator () (const Matrix<std::tuple<PixelT,
        PixelT, PixelT>> &m) const
    {
        uint vert_size = 2 * vert_radius + 1;
        uint hor_size = 2 * hor_radius + 1;
        ReturnT r = 0, g = 0, b = 0;
        KernelT sum_r = 0, sum_g = 0, sum_b = 0;
        for (uint i = 0; i < vert_size; ++i) {
            for (uint j = 0; j < hor_size; ++j) {
                // Tie is useful for taking elements from tuple
                std::tie(r, g, b) = m(i, j);
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
        return std::make_tuple(r, g, b);
    }
};

template<typename ReturnT, typename PixelT, typename KernelT>
Matrix<std::tuple<ReturnT, ReturnT, ReturnT>> custom(
    Matrix<std::tuple<PixelT, PixelT, PixelT>> &src_image, Matrix<KernelT> kernel, uint vr = 1,
    uint hr = 1) {
    return src_image.unary_map(KernelOp<ReturnT, PixelT, KernelT>(kernel, vr, hr));
}

template<typename PixelT>
Matrix<std::tuple<PixelT, PixelT, PixelT>> sobel_x(
    Matrix<std::tuple<uint, uint, uint>> &src_image) {
    Matrix<int> kernel = {{-1, 0, 1}};
    return custom<PixelT, uint, int>(src_image, kernel, 0, 1);
}

template<typename PixelT>
Matrix<std::tuple<PixelT, PixelT, PixelT>> sobel_y(
    Matrix<std::tuple<uint, uint, uint>> &src_image) {
    Matrix<int> kernel(3, 1);
    kernel(0, 0) = 1;
    kernel(1, 0) = 0;
    kernel(2, 0) = -1;
    return custom<PixelT, uint, int>(src_image, kernel, 1, 0);
}

Matrix<std::tuple<float, float>> edge_mapping(Matrix<std::tuple<uint, uint, uint>> &src_image) {
    Matrix<std::tuple<int, int, int>> sobelX = sobel_x<int>(src_image);
    Matrix<std::tuple<int, int, int>> sobelY = sobel_y<int>(src_image);
    Matrix<std::tuple<float, float>> result(src_image.n_rows, src_image.n_cols);
    for (uint row = 0; row < result.n_rows; ++row) {
        for (uint col = 0; col < result.n_cols; ++col) {
            int pixel_x = std::get<0>(sobelX(row, col));
            int pixel_y = std::get<0>(sobelY(row, col));
            float norm = std::sqrt(pixel_x * pixel_x + pixel_y * pixel_y);
            float angle = std::atan2(pixel_y, pixel_x);
            result(row, col) = std::make_tuple(norm, angle);
        }
    }
    return result;
}

// Load list of files and its labels from 'data_file' and
// stores it in 'file_list'
void LoadFileList(const string& data_file, TFileList* file_list) {
    ifstream stream(data_file.c_str());

    string filename;
    int label;
    
    int char_idx = data_file.size() - 1;
    for (; char_idx >= 0; --char_idx)
        if (data_file[char_idx] == '/' || data_file[char_idx] == '\\')
            break;
    string data_path = data_file.substr(0,char_idx+1);
    
    while(!stream.eof() && !stream.fail()) {
        stream >> filename >> label;
        if (filename.size())
            file_list->push_back(make_pair(data_path + filename, label));
    }

    stream.close();
}

// Load images by list of files 'file_list' and store them in 'data_set'
void LoadImages(const TFileList& file_list, TDataSet* data_set) {
    for (size_t img_idx = 0; img_idx < file_list.size(); ++img_idx) {
            // Create image
        BMP* image = new BMP();
            // Read image from file
        image->ReadFromFile(file_list[img_idx].first.c_str());
            // Add image and it's label to dataset
        data_set->push_back(make_pair(image, file_list[img_idx].second));
    }
}

// Save result of prediction to file
void SavePredictions(const TFileList& file_list,
                     const TLabels& labels, 
                     const string& prediction_file) {
        // Check that list of files and list of labels has equal size 
    assert(file_list.size() == labels.size());
        // Open 'prediction_file' for writing
    ofstream stream(prediction_file.c_str());

        // Write file names and labels to stream
    for (size_t image_idx = 0; image_idx < file_list.size(); ++image_idx)
        stream << file_list[image_idx].first << " " << labels[image_idx] << endl;
    stream.close();
}

void make_histogram(std::vector<float> *result, Matrix<std::tuple<float, float>> &region,
    uint sect_num) {
    result->resize(sect_num);
    for (uint row = 0; row < region.n_rows; ++row) {
        for (uint col = 0; col < region.n_cols; ++col) {
            float angle = std::get<1>(region(row, col));
            uint index = 0;
            if ((angle >= -(sect_num - 1) * PI / sect_num) ||
                (angle < (sect_num - 1) * PI / sect_num)) {
                index = std::min(std::max(static_cast<int>((angle + (sect_num - 1) *
                    PI / sect_num) / (2 * (sect_num - 1) * PI / sect_num) *
                    (sect_num - 1)) + 1, 1), static_cast<int>(sect_num) - 1);
            }
            (*result)[index] += std::get<0>(region(row, col));
        }
    }
    float sum = 0;
    for (auto vit = result->cbegin(); vit != result->cend(); ++vit) {
        sum += (*vit) * (*vit);
    }
    sum = std::sqrt(sum);
    if (sum > 0) {
        for (auto vit = result->begin(); vit != result->end(); ++vit) {
            *vit /= sum;
        }
    }
}

void LBP(std::vector<float> *hash, Matrix<std::tuple<uint, uint, uint>> &region) {
    hash->resize(256);
    uint row_size = region.n_rows - 1, col_size = region.n_cols - 1;
    for (uint row = 1; row < row_size; ++row) {
        for (uint col = 1; col < col_size; ++col) {
            uint value = 0, central_pixel = std::get<0>(region(row, col));
            for (uint i = 0; i < 3; ++i) {
                for (uint j = 0; j < 3; ++j) {
                    if ((i == 1) && (j == 1)) {
                        continue;
                    }
                    uint cmp = std::get<0>(region(row - 1 + i, col - 1 + j));
                    if (cmp >= central_pixel) {
                        value |= 1 << (3 * i + j);
                    }
                }
            }
            value = ((value & 480) >> 1) | (value & 15);
            (*hash)[value] += 1;
        }
    }
    float sum = 0;
    for (auto vit = hash->cbegin(); vit != hash->cend(); ++vit) {
        sum += (*vit) * (*vit);
    }
    sum = std::sqrt(sum);
    if (sum > 0) {
        for (auto vit = hash->begin(); vit != hash->cend(); ++vit) {
            *vit /= sum;
        }
    }
}

void colored_triple(std::vector<float> *mean_colors,
    Matrix<std::tuple<uint, uint, uint>> &region) {
    int mean_r = 0, mean_g = 0, mean_b = 0;
    for (uint row = 0; row < region.n_rows; ++row) {
        for (uint col = 0; col < region.n_cols; ++col) {
            int r = 0, g = 0, b = 0;
            std::tie(r, g, b) = region(row, col);
            mean_r += r;
            mean_g += g;
            mean_b += b;
        }
    }
    mean_colors->push_back(static_cast<float>(mean_r) / (255.0 * region.n_rows * region.n_cols));
    mean_colors->push_back(static_cast<float>(mean_g) / (255.0 * region.n_rows * region.n_cols));
    mean_colors->push_back(static_cast<float>(mean_b) / (255.0 * region.n_rows * region.n_cols));
}

struct HOG_args {
    uint vert_cell_num;
    uint hor_cell_num;
    uint sect_num;
    HOG_args(uint vcn = 16, uint hcn = 16, uint sn = 8): vert_cell_num(vcn), hor_cell_num(hcn),
    sect_num(sn) {}
};

void extract_HOG(Matrix<std::tuple<float, float>> &edge_map, vector<float> *one_image_features,
    HOG_args &args) {
    uint cell_height = edge_map.n_rows / args.vert_cell_num;
    uint cell_width = edge_map.n_cols / args.hor_cell_num;
    for (uint i = 0; i < args.vert_cell_num; ++i) {
        for (uint j = 0; j < args.hor_cell_num; ++j) {
            uint height = cell_height, width = cell_width;
            if (i == args.vert_cell_num - 1) {
                height = edge_map.n_rows - i * cell_height;
            }
            if (j == args.hor_cell_num - 1) {
                width = edge_map.n_cols - j * cell_width;
            }
            auto cell = edge_map.submatrix(i * cell_height, j * cell_width, height, width);
            std::vector<float> histogram;
            make_histogram(&histogram, cell, args.sect_num);
            for (uint vidx = 0; vidx < histogram.size(); ++vidx) {
                one_image_features->push_back(histogram[vidx]);
            }
        }
    }
}

struct LBP_args {
    uint LBP_vcn;
    uint LBP_hcn;
    LBP_args(uint vcn = 4, uint hcn = 4): LBP_vcn(vcn), LBP_hcn(hcn) {}
};

void extract_LBP(Matrix<std::tuple<uint, uint, uint>> &grayscale_image_extended,
    vector<float> *one_image_features, LBP_args &args) {
    uint LBP_ch = (grayscale_image_extended.n_rows - 2) / args.LBP_vcn;
    uint LBP_cw = (grayscale_image_extended.n_cols - 2) / args.LBP_hcn;
    for (uint i = 0; i < args.LBP_vcn; ++i) {
        for (uint j = 0; j < args.LBP_hcn; ++j) {
            uint height = LBP_ch + 2, width = LBP_cw + 2;
            if (i == args.LBP_vcn - 1) {
                height = grayscale_image_extended.n_rows - i * LBP_ch;
            }
            if (j == args.LBP_hcn - 1) {
                width = grayscale_image_extended.n_cols - j * LBP_cw;
            }
            auto hash_cell = grayscale_image_extended.submatrix(i * LBP_ch, j * LBP_cw, height,
                width);
            std::vector<float> hash;
            LBP(&hash, hash_cell);
            for (uint vidx = 0; vidx < hash.size(); ++vidx) {
                one_image_features->push_back(hash[vidx]);
            }
        }
    }
}

struct Color_args {
    uint color_vcn;
    uint color_hcn;
    Color_args(uint vcn = 4, uint hcn = 4): color_vcn(vcn), color_hcn(hcn) {}
};

void extract_Color(Matrix<std::tuple<uint, uint, uint>> &image, vector<float> *one_image_features,
    Color_args &args) {
    uint color_ch = image.n_rows / args.color_vcn;
    uint color_cw = image.n_cols / args.color_hcn;
    for (uint i = 0; i < args.color_vcn; ++i) {
        for (uint j = 0; j < args.color_hcn; ++j) {
            uint height = color_ch, width = color_cw;
            if (i == args.color_vcn - 1) {
                height = image.n_rows - i * color_ch;
            }
            if (j == args.color_hcn - 1) {
                width = image.n_cols - j * color_cw;
            }
            auto colored_cell = image.submatrix(i * color_ch, j * color_cw, height, width);
            std::vector<float> mean_pixel;
            colored_triple(&mean_pixel, colored_cell);
            for (uint vidx = 0; vidx < mean_pixel.size(); ++vidx) {
                one_image_features->push_back(mean_pixel[vidx]);
            }
        }
    }
}

// Exatract features from dataset.
// You should implement this function by yourself =)
void ExtractFeatures(const TDataSet& data_set, TFeatures* features) {
    HOG_args hog_arg(8, 8, 16);
    LBP_args lbp_arg(4, 4);
    Color_args color_arg(8, 8);
    for (size_t image_idx = 0; image_idx < data_set.size(); ++image_idx) {
        
        // PLACE YOUR CODE HERE
        Matrix<std::tuple<uint, uint, uint>> image = fromBMPtoMatrix(data_set[image_idx].first);
        Matrix<std::tuple<uint, uint, uint>> grayscale_image = grayscale(image);
        Matrix<std::tuple<uint, uint, uint>> grayscale_image_extended =
            grayscale_image.extra_borders(1, 1);
        Matrix<std::tuple<float, float>> edge_map = edge_mapping(grayscale_image);
        vector<float> one_image_features;
        extract_HOG(edge_map, &one_image_features, hog_arg);
        extract_LBP(grayscale_image_extended, &one_image_features, lbp_arg);
        extract_Color(image, &one_image_features, color_arg);
        features->push_back(std::make_pair(one_image_features, data_set[image_idx].second));
    }
}

// Clear dataset structure
void ClearDataset(TDataSet* data_set) {
        // Delete all images from dataset
    for (size_t image_idx = 0; image_idx < data_set->size(); ++image_idx)
        delete (*data_set)[image_idx].first;
        // Clear dataset
    data_set->clear();
}

// Train SVM classifier using data from 'data_file' and save trained model
// to 'model_file'
void TrainClassifier(const string& data_file, const string& model_file) {
        // List of image file names and its labels
    TFileList file_list;
        // Structure of images and its labels
    TDataSet data_set;
        // Structure of features of images and its labels
    TFeatures features;
        // Model which would be trained
    TModel model;
        // Parameters of classifier
    TClassifierParams params;
    
        // Load list of image file names and its labels
    LoadFileList(data_file, &file_list);
        // Load images
    LoadImages(file_list, &data_set);
        // Extract features from images
    ExtractFeatures(data_set, &features);

        // PLACE YOUR CODE HERE
        // You can change parameters of classifier here
    params.C = 0.01;
    TClassifier classifier(params);
        // Train classifier
    classifier.Train(features, &model);
        // Save model to file
    model.Save(model_file);
        // Clear dataset structure
    ClearDataset(&data_set);
}

// Predict data from 'data_file' using model from 'model_file' and
// save predictions to 'prediction_file'
void PredictData(const string& data_file,
                 const string& model_file,
                 const string& prediction_file) {
        // List of image file names and its labels
    TFileList file_list;
        // Structure of images and its labels
    TDataSet data_set;
        // Structure of features of images and its labels
    TFeatures features;
        // List of image labels
    TLabels labels;

        // Load list of image file names and its labels
    LoadFileList(data_file, &file_list);
        // Load images
    LoadImages(file_list, &data_set);
        // Extract features from images
    ExtractFeatures(data_set, &features);

        // Classifier 
    TClassifier classifier = TClassifier(TClassifierParams());
        // Trained model
    TModel model;
        // Load model from file
    model.Load(model_file);
        // Predict images by its features using 'model' and store predictions
        // to 'labels'
    classifier.Predict(features, model, &labels);

        // Save predictions
    SavePredictions(file_list, labels, prediction_file);
        // Clear dataset structure
    ClearDataset(&data_set);
}

int main(int argc, char** argv) {
    // Command line options parser
    ArgvParser cmd;
        // Description of program
    cmd.setIntroductoryDescription("Machine graphics course, task 2. CMC MSU, 2014.");
        // Add help option
    cmd.setHelpOption("h", "help", "Print this help message");
        // Add other options
    cmd.defineOption("data_set", "File with dataset",
        ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
    cmd.defineOption("model", "Path to file to save or load model",
        ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
    cmd.defineOption("predicted_labels", "Path to file to save prediction results",
        ArgvParser::OptionRequiresValue);
    cmd.defineOption("train", "Train classifier");
    cmd.defineOption("predict", "Predict dataset");
        
        // Add options aliases
    cmd.defineOptionAlternative("data_set", "d");
    cmd.defineOptionAlternative("model", "m");
    cmd.defineOptionAlternative("predicted_labels", "l");
    cmd.defineOptionAlternative("train", "t");
    cmd.defineOptionAlternative("predict", "p");

        // Parse options
    int result = cmd.parse(argc, argv);

        // Check for errors or help option
    if (result) {
        cout << cmd.parseErrorDescription(result) << endl;
        return result;
    }

        // Get values 
    string data_file = cmd.optionValue("data_set");
    string model_file = cmd.optionValue("model");
    bool train = cmd.foundOption("train");
    bool predict = cmd.foundOption("predict");

        // If we need to train classifier
    if (train)
        TrainClassifier(data_file, model_file);
        // If we need to predict data
    if (predict) {
            // You must declare file to save images
        if (!cmd.foundOption("predicted_labels")) {
            cerr << "Error! Option --predicted_labels not found!" << endl;
            return 1;
        }
            // File to save predictions
        string prediction_file = cmd.optionValue("predicted_labels");
            // Predict data
        PredictData(data_file, model_file, prediction_file);
    }
}
