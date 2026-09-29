#include <iostream>

#include "line.h"

int main() {
    Webbing webbing{5.0, 10.0, 50.0};
    Segment segment(webbing, webbing, 10.0, 10.5);

    const int num_segments = 10;
    Line line(num_segments, std::vector<Segment>(num_segments, segment));
    std::cout << "Line with " << line.num_segments() << " segments\n";
    return 0;
}
