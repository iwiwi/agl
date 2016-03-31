#pragma once
#include "two_edge_cc_filter.h"
#include "cut_tree_with_2ecc.h"
#include "dinitz.h"
#include "bi_dinitz.h"
#include "plain_gusfield.h"

namespace agl {
using cut_tree = agl::cut_tree_internal::two_edge_cc_filter<cut_tree_with_2ecc>;
} //namespace agl