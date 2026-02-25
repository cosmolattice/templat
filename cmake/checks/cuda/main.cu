__global__ void hello() {}
int main() { hello<<<1,1>>>(); return 0; }
