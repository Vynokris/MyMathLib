#pragma once
#include "Vector2.h"
#include "Vector3.h"
constexpr rsize_t MAX_BONE_WEIGHTS = 4;

namespace Maths
{
    // - Vertex: Point in 3D space with rendering data - //
    struct Vertex
    {
        Vector3 pos;    // Vertex position.
        Vector2 uv;     // Vertex texture coordinates.
        Vector3 normal; // Vertex normal vector.

        bool operator==(const Vertex& other) const
        {
            return pos    == other.pos
                && uv     == other.uv
                && normal == other.normal;
        }
    };

    // - TangentVertex: Point in 3D space with rendering data and tangent space info - //
    struct TangentVertex
    {
        Vector3 pos;       // Vertex position.
        Vector2 uv;        // Vertex texture coordinates.
        Vector3 normal;    // Vertex normal vector.
        Vector3 tangent;   // Vertex tangent vector.
        Vector3 bitangent; // Vertex bitangent vector (orthogonal to normal and tangent).

        bool operator==(const TangentVertex& other) const
        {
            return pos       == other.pos
                && uv        == other.uv
                && normal    == other.normal
                && tangent   == other.tangent
                && bitangent == other.bitangent;
        }
    };

    struct AnimatedVertex : TangentVertex
    {
        int   boneIDs[MAX_BONE_WEIGHTS] = {}; // Bone indices which influence this vertex.
        float weights[MAX_BONE_WEIGHTS] = {}; // Weights that dictate how much each bone influences this vertex.

        AnimatedVertex()
        {
            for (int i = 0; i < 4; i++)
            {
                boneIDs[i] = -1;
                weights[i] = 0.0f;
            }
        }
    };

    // - VertexIndices: Holds indices for vertex data - //
    struct VertexIndices
    {
        uint32_t pos;    // Position index.
        uint32_t uv;     // Texture coordinate index.
        uint32_t normal; // Normal index.

        bool operator==(const VertexIndices& other) const
        {
            return pos    == other.pos
                && uv     == other.uv
                && normal == other.normal;
        }
    };
}
