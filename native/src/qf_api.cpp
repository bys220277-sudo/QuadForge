#include "qf_api.h"

#include "qf_pipeline.h"

using namespace qf;

namespace {
struct Handle {
    Result res;
    std::string log;
};
}  // namespace

QF_API int qf_api_version(void) { return QF_API_VERSION; }

QF_API int qf_default_params(QFParams* p) {
    if (!p) return -1;
    std::memset(p, 0, sizeof(*p));
    p->struct_size = (int32_t)sizeof(QFParams);
    p->target_faces = 5000;
    p->sharp_angle = 0.6108652381980153;  // 35°
    p->detect_sharp = 1;
    p->preserve_boundary = 1;
    p->adaptive = 0.0;
    p->curvature_align = 0.5;
    p->guide_strength = 0.7;
    p->density_strength = 0.5;
    p->seed = 1;
    p->quality = 2;
    p->relax_iterations = 0;
    return 0;
}

QF_API void qf_set_threads(int n) { set_num_threads(n); }
QF_API int qf_get_threads(void) { return num_threads(); }

QF_API void* qf_run(const QFParams* p, const QFInput* in, QFStatus* st) {
    Handle* h = new Handle();
    Ctx ctx;
    ctx.st = reinterpret_cast<Status*>(st);
    ctx.verbose = p && p->verbose;
    try {
        if (!p || p->struct_size != (int32_t)sizeof(QFParams)) {
            h->res.error = "несовпадение версии структуры параметров";
        } else if (!in || in->struct_size != (int32_t)sizeof(QFInput)) {
            h->res.error = "несовпадение версии структуры входа";
        } else {
            if (p->threads > 0 && p->threads != num_threads()) set_num_threads(p->threads);
            run_pipeline(*p, *in, ctx, h->res);
        }
    } catch (const Cancelled&) {
        h->res.ok = false;
        h->res.error = "отменено";
    } catch (const std::exception& e) {
        h->res.ok = false;
        h->res.error = std::string("исключение: ") + e.what();
    } catch (...) {
        h->res.ok = false;
        h->res.error = "неизвестное исключение";
    }
    h->log = ctx.log;
    if (st) st->done = 1;
    return h;
}

QF_API int qf_result_ok(void* r) { return r && static_cast<Handle*>(r)->res.ok ? 1 : 0; }
QF_API const char* qf_result_error(void* r) { return r ? static_cast<Handle*>(r)->res.error.c_str() : ""; }
QF_API const char* qf_result_log(void* r) { return r ? static_cast<Handle*>(r)->log.c_str() : ""; }

QF_API int qf_result_mesh(void* r, const float** verts, int* nverts, const int** face_sizes,
                          int* nfaces, const int** indices, int* nindices) {
    if (!r) return -1;
    Result& R = static_cast<Handle*>(r)->res;
    *verts = R.V.data();
    *nverts = (int)(R.V.size() / 3);
    *face_sizes = R.face_sizes.data();
    *nfaces = (int)R.face_sizes.size();
    *indices = R.indices.data();
    *nindices = (int)R.indices.size();
    return 0;
}

QF_API int qf_result_array(void* r, const char* name, const void** data, int64_t* n, int* dtype) {
    if (!r || !name) return -1;
    Result& R = static_cast<Handle*>(r)->res;
    auto it = R.darr.find(name);
    if (it != R.darr.end()) {
        *data = it->second.data();
        *n = (int64_t)it->second.size();
        *dtype = 0;
        return 0;
    }
    auto it2 = R.iarr.find(name);
    if (it2 != R.iarr.end()) {
        *data = it2->second.data();
        *n = (int64_t)it2->second.size();
        *dtype = 1;
        return 0;
    }
    return -1;
}

QF_API void qf_result_free(void* r) { delete static_cast<Handle*>(r); }
