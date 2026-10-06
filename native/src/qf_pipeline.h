#pragma once

#include "qf_api.h"
#include "qf_common.h"

#include <map>

namespace qf {

struct Result {
    bool ok = false;
    std::string error;
    std::vector<float> V;
    std::vector<int> face_sizes;
    std::vector<int> indices;
    std::map<std::string, std::vector<double>> darr;
    std::map<std::string, std::vector<int>> iarr;
};

void run_pipeline(const QFParams& p, const QFInput& in, Ctx& ctx, Result& res);

}  // namespace qf
