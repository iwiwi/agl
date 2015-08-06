#pragma once

#include "radix_heap.h"
#include "dijkstra_heap.h"

namespace agl {
template<typename KeyT>
using heap = radix_heap::radix_heap<KeyT>;

template<typename KeyT, typename ValueT>
using pair_heap = radix_heap::pair_radix_heap<KeyT, ValueT>;
}
