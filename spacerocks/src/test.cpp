// Your First C++ Program

#include <iostream>

int main() {

  // double a [1000000000] = { };
  for (int i = 0; i < 10; i++) {
    // a[i] = i / (2.0 / 100.0 - (0.1 * 0.1 + 0.1 * 0.1 + 0.1 * 0.1) / 1.0);
    std::cout << rand() / double(RAND_MAX) << "\n";
    }
    return 0;
}
