/*switch (channelA) {
    case 2: {//srcA = B, srcB = G
        for (int row = 0; row < find_pos.rows; ++row) {
            for (int col = 0; col < find_pos.cols; ++col) {
                std::get<2>(result(row, col)) = std::get<2>(srcA(row + A_row_shift, col + A_col_shift));
                std::get<1>(result(row, col)) = std::get<1>(srcB(row - B_row_shift, col - B_col_shift));
                std::get<0>(result(row, col)) = std::get<0>(srcA(row + A_row_shift, col + A_col_shift));
            }
        }
        break;
    }
    case 0: {//srcA = R, srcB = BG
        for (int row = 0; row < find_pos.rows; ++row) {
            for (int col = 0; col < find_pos.cols; ++col) {
                std::get<2>(result(row, col)) = std::get<2>(srcB(row - B_row_shift, col - B_col_shift));
                std::get<1>(result(row, col)) = std::get<1>(srcB(row - B_row_shift, col - B_col_shift));
                std::get<0>(result(row, col)) = std::get<0>(srcA(row + A_row_shift, col + A_col_shift));
            }
        }
        break;
    }
    default: break;
}*/
/*Image BG = mse(channels.B, channels.G, delta_shift, 2);
Image BGR = mse(channels.R, BG, delta_shift, 0);*/
Image mse(Image &srcA, Image &srcB, int delta_shift, uint channelA)
{
    cout << "entered mse\n";
    Pos find_pos;
    Image Big;
    Image Little; 
    if (srcA.n_rows < srcB.n_rows) {
        find_pos.first_big = false;
        Big = srcB;
        Little = srcA;
    } else {
        Big = srcA;
        Little = srcB;
    }
    double min_error = 0.0;
    for (int i = -delta_shift; i < delta_shift + 1; ++i) {
        int rows = Little.n_rows;
        if (i < 0) {
            rows += i;
        } else {
            if (i > static_cast<int>(Big.n_rows - Little.n_rows)) {
                rows += Big.n_rows - Little.n_rows - i;
            }
        }
        for (int j = -delta_shift; j < delta_shift + 1; ++j) {
            double error = 0.0;
            bool div_flag = true;
            /*int rows = Little.n_rows;
            if (i < 0) {
                rows += i;
            } else {
                if (i > static_cast<int>(Big.n_rows - Little.n_rows)) {
                    rows += Big.n_rows - Little.n_rows - i;
                }
            }*/
            int cols = Little.n_cols;
            if (j < 0) {
                cols += j;
            } else {
                cols -= j;
            }
            int big_col_shift = ((j > 0) ? j : 0), little_col_shift = ((j < 0) ? j : 0);
            for (int row = 0; row < rows; ++row) {
                for (int col = 0; col < cols; ++col) {
                    error += std::pow(std::get<1>(Big(row + ((i > 0) ? i : 0), col + big_col_shift)) - std::get<1>(Little(row - ((i < 0) ? i : 0), col - little_col_shift)), 2);
                    if (div_flag && (error > rows * cols)) {
                        div_flag = false;
                        error /= rows * cols;
                    }
                }
            }
            if (div_flag) {
                error /= rows * cols;
            }
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
    int big_row_shift = ((find_pos.row_shift > 0) ? find_pos.row_shift : 0), little_row_shift = ((find_pos.row_shift < 0) ? find_pos.row_shift : 0);
    int big_col_shift = ((find_pos.col_shift > 0) ? find_pos.col_shift : 0), little_col_shift = ((find_pos.col_shift < 0) ? find_pos.col_shift : 0);
    if (find_pos.first_big) {
        switch (channelA) {
            case 2: {//Big = B(1), Little = G(1)
                for (int row = 0; row < find_pos.rows; ++row) {
                    for (int col = 0; col < find_pos.cols; ++col) {
                        std::get<2>(result(row, col)) = std::get<2>(Big(row + big_row_shift, col + big_col_shift));
                        std::get<1>(result(row, col)) = std::get<1>(Little(row - little_row_shift, col - little_col_shift));
                    }
                }
                break;
            }
            case 0: {//Big = R(1), Little = BG(1)
                for (int row = 0; row < find_pos.rows; ++row) {
                    for (int col = 0; col < find_pos.cols; ++col) {
                        std::get<2>(result(row, col)) = std::get<2>(Little(row - little_row_shift, col - little_col_shift));
                        std::get<1>(result(row, col)) = std::get<1>(Little(row - little_row_shift, col - little_col_shift));
                        std::get<0>(result(row, col)) = std::get<0>(Big(row + big_row_shift, col + big_col_shift));
                    }
                }
                break;
            }
            default: break;
        }
    } else {
        switch (channelA) {
            case 2: {//Big = G(1), Little = B(1)
                for (int row = 0; row < find_pos.rows; ++row) {
                    for (int col = 0; col < find_pos.cols; ++col) {
                        std::get<2>(result(row, col)) = std::get<2>(Little(row - little_row_shift, col - little_col_shift));
                        std::get<1>(result(row, col)) = std::get<1>(Big(row + big_row_shift, col + big_col_shift));
                    }
                }
                break;
            }
            case 0: {//Big = BG(1), Little = R(1)
                for (int row = 0; row < find_pos.rows; ++row) {
                    for (int col = 0; col < find_pos.cols; ++col) {
                        std::get<2>(result(row, col)) = std::get<2>(Big(row + big_row_shift, col + big_col_shift));
                        std::get<1>(result(row, col)) = std::get<1>(Big(row + big_row_shift, col + big_col_shift));
                        std::get<0>(result(row, col)) = std::get<0>(Little(row - little_row_shift, col - little_col_shift));
                    }
                }
                break;
            }
            default: break;
        }
    }
    return result;
}
/*
    //type 0
    for (uint i = 0; i < errors.size(); ++i) {
        for (uint j = 0; j < errors[0].size(); ++j) {
            for (uint row_i = 0; row_i < new_rows - i; ++row_i) {
                for (uint col_j = 0; col_j < new_cols - j; ++col_j) {
                    error += (std::get<1>(srcA(i + row_i, j + col_j)) - std::get<1>(srcB(row_i, col_j))) * (std::get<1>(srcA(i + row_i, j + col_j)) - std::get<1>(srcB(row_i, col_j)));
                    if (div_flag && error > new_rows * new_cols) {
                        div_flag = false;
                        error /= new_rows * new_cols;
                    }
                }
            }
            errors[i][j] = error;
            error = 0.0;
        }
    }
    min_error = errors[0][0];
    for (uint i = 0; i < errors.size(); ++i) {
        for (uint j = 0; j < errors[0].size(); ++j) {
            if (errors[i][j] < min_error) {
                min_error = errors[i][j];
                find_pos.row = i;
                find_pos.col = j;
            }
        }
    }
    //type 1
    div_flag = true;
    for (uint i = 0; i < errors.size(); ++i) {
        for (uint j = 0; j < errors[0].size(); ++j) {
            for (uint row_i = 0; row_i < new_rows - i; ++row_i) {
                for (uint col_j = 0; col_j < new_cols - j; ++col_j) {
                    error += (std::get<1>(srcB(i + row_i, j + col_j)) - std::get<1>(srcA(row_i, col_j))) * (std::get<1>(srcB(i + row_i, j + col_j)) - std::get<1>(srcA(row_i, col_j)));
                    if (div_flag && error > new_rows * new_cols) {
                        div_flag = false;
                        error /= new_rows * new_cols;
                    }
                }
            }
            errors[i][j] = error;
            error = 0.0;
        }
    }
    for (uint i = 0; i < errors.size(); ++i) {
        for (uint j = 0; j < errors[0].size(); ++j) {
            if (errors[i][j] < min_error) {
                min_error = errors[i][j];
                find_pos.row = i;
                find_pos.col = j;
                find_pos.type = 1;
            }
        }
    }
    //type 2
    div_flag = true;
    for (uint i = 0; i < errors.size(); ++i) {
        for (uint j = 0; j < errors[0].size(); ++j) {
            for (uint row_i = 0; row_i < new_rows - i; ++row_i) {
                for (uint col_j = 0; col_j < new_cols - j; ++col_j) {
                    error += (std::get<1>(srcA(i + row_i, col_j)) - std::get<1>(srcB(row_i, j + col_j))) * (std::get<1>(srcA(i + row_i, col_j)) - std::get<1>(srcB(row_i, j + col_j)));
                    if (div_flag && error > new_rows * new_cols) {
                        div_flag = false;
                        error /= new_rows * new_cols;
                    }
                }
            }
            errors[i][j] = error;
            error = 0.0;
        }
    }
    for (uint i = 0; i < errors.size(); ++i) {
        for (uint j = 0; j < errors[0].size(); ++j) {
            if (errors[i][j] < min_error) {
                min_error = errors[i][j];
                find_pos.row = i;
                find_pos.col = j;
                find_pos.type = 2;
            }
        }
    }
    //type 3
    div_flag = true;
    for (uint i = 0; i < errors.size(); ++i) {
        for (uint j = 0; j < errors[0].size(); ++j) {
            for (uint row_i = 0; row_i < new_rows - i; ++row_i) {
                for (uint col_j = 0; col_j < new_cols - j; ++col_j) {
                    error += (std::get<1>(srcB(i + row_i, col_j)) - std::get<1>(srcA(row_i, j + col_j))) * (std::get<1>(srcB(i + row_i, col_j)) - std::get<1>(srcA(row_i, j + col_j)));
                    if (div_flag && error > new_rows * new_cols) {
                        div_flag = false;
                        error /= new_rows * new_cols;
                    }
                }
            }
            errors[i][j] = error;
            error = 0.0;
        }
    }
    for (uint i = 0; i < errors.size(); ++i) {
        for (uint j = 0; j < errors[0].size(); ++j) {
            if (errors[i][j] < min_error) {
                min_error = errors[i][j];
                find_pos.row = i;
                find_pos.col = j;
                find_pos.type = 3;
            }
        }
    }
    //find the best shift
    Image result(new_rows - find_pos.row, new_cols - find_pos.col);
    if (cA == 2) {
        switch (find_pos.type) {
            case 0: {
                for (uint row = 0; row < result.n_rows; ++row) {
                    for (uint col = 0; col < result.n_cols; ++col) {
                        std::get<2>(result(row, col)) = std::get<2>(srcA(find_pos.row + row, find_pos.col + col));
                        std::get<1>(result(row, col)) = std::get<1>(srcB(row, col));
                    }
                }
                break;
            }
            case 1: {
                for (uint row = 0; row < result.n_rows; ++row) {
                    for (uint col = 0; col < result.n_cols; ++col) {
                        std::get<1>(result(row, col)) = std::get<1>(srcB(find_pos.row + row, find_pos.col + col));
                        std::get<2>(result(row, col)) = std::get<2>(srcA(row, col));
                    }
                }
                break;
            }
            case 2: {
                for (uint row = 0; row < result.n_rows; ++row) {
                    for (uint col = 0; col < result.n_cols; ++col) {
                        std::get<2>(result(row, col)) = std::get<2>(srcA(find_pos.row + row, col));
                        std::get<1>(result(row, col)) = std::get<1>(srcB(row, find_pos.col + col));
                    }
                }
                break;
            }
            case 3: {
                for (uint row = 0; row < result.n_rows; ++row) {
                    for (uint col = 0; col < result.n_cols; ++col) {
                        std::get<1>(result(row, col)) = std::get<1>(srcB(find_pos.row + row, col));
                        std::get<2>(result(row, col)) = std::get<2>(srcA(row, find_pos.col + col));
                    }
                }
                break;
            }
            default: 
                break;
        }
    } else {
        switch (find_pos.type) {
            case 0: {
                for (uint row = 0; row < result.n_rows; ++row) {
                    for (uint col = 0; col < result.n_cols; ++col) {
                        std::get<2>(result(row, col)) = std::get<2>(srcA(find_pos.row + row, find_pos.col + col));
                        std::get<1>(result(row, col)) = std::get<1>(srcA(find_pos.row + row, find_pos.col + col));
                        std::get<0>(result(row, col)) = std::get<0>(srcB(row, col));
                    }
                }
                break;
            }
            case 1: {
                for (uint row = 0; row < result.n_rows; ++row) {
                    for (uint col = 0; col < result.n_cols; ++col) {
                        std::get<2>(result(row, col)) = std::get<2>(srcA(row, col));
                        std::get<0>(result(row, col)) = std::get<0>(srcB(find_pos.row + row, find_pos.col + col));
                        std::get<1>(result(row, col)) = std::get<1>(srcA(row, col));
                    }
                }
                break;
            }
            case 2: {
                for (uint row = 0; row < result.n_rows; ++row) {
                    for (uint col = 0; col < result.n_cols; ++col) {
                        std::get<2>(result(row, col)) = std::get<2>(srcA(find_pos.row + row, col));
                        std::get<1>(result(row, col)) = std::get<1>(srcA(find_pos.row + row, col));
                        std::get<0>(result(row, col)) = std::get<0>(srcB(row, find_pos.col + col));
                    }
                }
                break;
            }
            case 3: {
                for (uint row = 0; row < result.n_rows; ++row) {
                    for (uint col = 0; col < result.n_cols; ++col) {
                        std::get<1>(result(row, col)) = std::get<1>(srcA(row, find_pos.col + col));
                        std::get<0>(result(row, col)) = std::get<0>(srcB(find_pos.row + row, col));
                        std::get<1>(result(row, col)) = std::get<1>(srcA(row, find_pos.col + col));
                    }
                }
                break;
            }
            default: 
                break;
        }
    }
    */