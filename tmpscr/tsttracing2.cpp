#include <cstdio>
#include <iostream>

#define TRACED(name, ...) name(__VA_ARGS__) { \
      std::printf("Calling %s\n", #name);


int TRACED(add, int a, int b)
   return a + b;
}

static void fooBar() {
    std::cout << "Calling " << __func__;
}

int main() {
   int res = add(3, 4);
   fooBar();
   return 0;
}
