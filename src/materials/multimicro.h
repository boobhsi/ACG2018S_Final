#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MATERIALS_MULTIMICRO_H
#define PBRT_MATERIALS_MULTIMICRO_H

#include "pbrt.h"
#include "material.h"
#include "acg_final/MicrosurfaceScattering.h"
#include <memory>

namespace pbrt {

// PlasticMaterial Declarations
class MultiMicroMaterial : public Material {
  public:
    // PlasticMaterial Public Methods
    MultiMicroMaterial(const std::shared_ptr<Texture<Float>> &roughnessX,
                    const std::shared_ptr<Texture<Float>> &roughnessY,
                    const std::shared_ptr<Texture<Float>> &bumpMap,
                    bool uni, bool beck) :
          roughnessX(roughnessX),
          roughnessY(roughnessY),
          bumpMap(bumpMap),
          uni(uni),
          beck(beck) {
            std::shared_ptr<Microsurface> temp(new MicrosurfaceDielectric(uni, beck));
            m_microsurface = temp;
          }
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    std::shared_ptr<Texture<Float>> roughnessX, roughnessY, bumpMap;
    std::shared_ptr<Microsurface> m_microsurface;
    bool uni, beck;
};

MultiMicroMaterial *CreateMultiMicroMaterial(const TextureParams &mp);

}  // namespace pbrt

#endif