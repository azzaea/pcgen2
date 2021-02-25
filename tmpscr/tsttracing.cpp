#include <iostream>

static void fooBar() {
    std::cout << "Calling " << __func__;
}

int main() {
  fooBar();
}
