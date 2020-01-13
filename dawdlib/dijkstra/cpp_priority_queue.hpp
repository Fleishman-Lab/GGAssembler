#include <queue>
template <class T>
using greater_priority_queue = std::priority_queue<T, std::vector<T>, std::greater<T>>;