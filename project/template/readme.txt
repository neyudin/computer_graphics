Реализованные функции:
    Matrix<std::tuple<uint, uint, uint>> fromBMPtoMatrix(BMP *src_image);

    Matrix<std::tuple<uint, uint, uint>> grayscale(Matrix<std::tuple<uint, uint, uint>> &src_image);

    template<typename ReturnT, typename PixelT, typename KernelT>
        Matrix<std::tuple<ReturnT, ReturnT, ReturnT>> custom(
        Matrix<std::tuple<PixelT, PixelT, PixelT>> &src_image, Matrix<KernelT> kernel, uint vr = 1,
        uint hr = 1);

    template<typename PixelT>
        Matrix<std::tuple<PixelT, PixelT, PixelT>> sobel_x(
        Matrix<std::tuple<uint, uint, uint>> &src_image);

    template<typename PixelT>
        Matrix<std::tuple<PixelT, PixelT, PixelT>> sobel_y(
        Matrix<std::tuple<uint, uint, uint>> &src_image);

    Matrix<std::tuple<float, float>> edge_mapping(Matrix<std::tuple<uint, uint, uint>> &src_image);

    void make_histogram(std::vector<float> *result, Matrix<std::tuple<float, float>> &region,
    uint sect_num);

    void LBP(std::vector<float> *hash, Matrix<std::tuple<uint, uint, uint>> &region);

    void colored_triple(std::vector<float> *mean_colors,
        Matrix<std::tuple<uint, uint, uint>> &region);

    void extract_HOG(Matrix<std::tuple<float, float>> *edge_map, vector<float> *one_image_features,
        HOG_args *args);

    void extract_LBP(Matrix<std::tuple<uint, uint, uint>> *grayscale_image_extended,
        vector<float> *one_image_features, LBP_args *args);

    void extract_Color(Matrix<std::tuple<uint, uint, uint>> *image, vector<float> *one_image_features,
        Color_args *args);

    void ExtractFeatures(const TDataSet& data_set, TFeatures* features);
Реализованные структуры данных:
    template<typename ReturnT, typename PixelT, typename KernelT> class KernelOp;

    struct HOG_args;

    struct LBP_args;

    struct Color_args;
В ходе выполнения работы была сделана базовая часть и выполнены дополнительные задания: Локальные бинарные шаблоны, Цветовые признаки.
Мною были исправлены ошибки в каркасе, связанные с утечкой памяти, детектированной в valgrind, основываясь на советах на
https://courses.graphics.cs.msu.ru/mod/forum/discuss.php?d=39 . Была добавлена свместимость с компилятором g++-4.8 соответствующей командой в
Makefale.
Файл с обученной моделью: model.txt.
Выполненные дополнительные части задания:
1) Локальные бинарные шаблоны
    Реализация в функции: extract_LBP.
    Количество клеток: 16 (деление изображения на одинаковое число сегментов по горизонтали и по вертикали)
    Локальные бинарные шаблоны создаются по изображению в оттенках серого.
    Для построения дескриптора локальных бинарных шаблонов (Local Binary Patterns, LBP) исходное изображение было зеркально расширено
    по границам с помощью метода extra_borders(1, 1) класса Matrix<std::tuple<uint, uint, uint>>. Это было сделано для того, чтобы в
    клетках, на которые исходное изображение делится, каждый пиксель имел соседа. Вычисляются высота и ширина каждой клетки относительно
    исходного изображения. Эти значения используются для извлечения клеток из расширенного изображения так, чтобы каждая клетка раширенного
    изображения была клеткой исходного изображения с обрамляющими её границами шириной 1 пиксель. Каждый пиксель клетки сравнивается cо своими
    соседями. Соседи для каждого пикселя перечисляются одним и тем же образом: в двойном цикле, внешний - по строкам, внутренний - по столбцам.
    Получается вектор из 8 бинарных значений: 0 означает, что значение пикселя больше, чем значение соседа, а 1 - что не больше. Полученный
    вектор трактуется как бинарная запись числа от 0 до 255. Затем, для каждой клетки строится гистограмма, показывающая сколько раз каждое
    такое "число" встречается в этой клетке. Получится вектор из 256 значений, который нормализуется. Каждый такой вектор строится функцией
    LBP. Все гистограммы конкатенируются. Получившийся дескриптор локальных бинарных шаблонов конкатенируется с дескриптором HOG.
2) Цветовые признаки
    Реализация в функции: extract_Color.
    Количество клеток: 64 (деление изображения на одинаковое число сегментов по горизонтали и по вертикали)
    Цветовые признаки создаются по цветному изображению.
    Для построения дескриптора цветовых признаков исходное изображение делится на клетки. В каждой клетке считается средний цвет: три значения -
    среднее для каждого из каналов. Все значения приводятся к отрезку [0, 1] (делятся на 255) и вытягиваются в один вектор. Приведённые к
    отрезку [0, 1] средние значения для каждого из каналов в клетке вычисляются функцией colored_triple. Полученный вектор конкатенируется с
    дескриптором локальных бинарных шаблонов и дескриптором HOG.
