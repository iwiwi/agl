#pragma once
#include "two_edge_cc_filter.h"
#include "cut_tree_with_2ecc.h"
#include "dinitz.h"
#include "bi_dinitz.h"
#include "plain_gusfield.h"

namespace agl {
using cut_tree = agl::cut_tree_internal::two_edge_cc_filter<cut_tree_with_2ecc>; // fastest

using plain_gusfield_bi_dinitz = agl::cut_tree_internal::plain_gusfield<bi_dinitz>; //faster than plain_gusfield_dinitz
using plain_gusfield_dinitz = agl::cut_tree_internal::plain_gusfield<agl::cut_tree_internal::dinitz>;
} //namespace agl
