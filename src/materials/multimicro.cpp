
#include "materials/multimicro.h"
#include "interaction.h"
#include "reflection.h"
#include "texture.h"
#include "paramset.h"
#include <memory>

namespace pbrt {

void MultiMicroMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
                                                    MemoryArena &arena,
                                                    TransportMode mode,
                                                    bool allowMultipleLobes) const {
    if(bumpMap) Bump(bumpMap, si);
    float eta = 1.5f;
    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, eta);
    Float roughu = roughnessX->Evaluate(*si);
    Float roughv = roughnessY->Evaluate(*si);
    si->bsdf->Add(ARENA_ALLOC(arena, MultiMicroBSDF)(m_microsurface, roughu, roughv));
}

MultiMicroMaterial *CreateMultiMicroMaterial(const TextureParams &mp) {
    std::shared_ptr<Texture<Float>> roughnessX =
        mp.GetFloatTexture("roughnessX", .1f);
    std::shared_ptr<Texture<Float>> roughnessY =
        mp.GetFloatTexture("roughnessY", .1f);
    std::shared_ptr<Texture<Float>> bumpMap =
        mp.GetFloatTextureOrNull("bumpmap");
    bool uni = mp.FindBool("uniform", false);
    bool beck = mp.FindBool("beckmann", false);
    return new MultiMicroMaterial(roughnessX, roughnessY, bumpMap, uni, beck);
}

}  // namespace pbrt
