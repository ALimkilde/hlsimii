#include "physics.h"
 
LineModel::LineModel(const DiscreteLine& line, Vec2 anchor_left, Vec2 anchor_right, const Params& p)
   : line_(std::move(line)),
     anchor_left_(anchor_left),
     anchor_right_(anchor_right),
     params_(p),
     Fe_(2, line.num_elements())
{
}
