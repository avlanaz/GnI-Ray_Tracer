using System;

namespace RayTracer
{
    /// <summary>
    /// Class to represent a triangle in a scene represented by three vertices.
    /// </summary>
    public class Triangle : SceneEntity
    {
        private Vector3 v0, v1, v2;
        private Material material;

        /// <summary>
        /// Construct a triangle object given three vertices.
        /// </summary>
        /// <param name="v0">First vertex position</param>
        /// <param name="v1">Second vertex position</param>
        /// <param name="v2">Third vertex position</param>
        /// <param name="material">Material assigned to the triangle</param>
        public Triangle(Vector3 v0, Vector3 v1, Vector3 v2, Material material)
        {
            this.v0 = v0;
            this.v1 = v1;
            this.v2 = v2;
            this.material = material;
        }

        /// <summary>
        /// Determine if a ray intersects with the triangle, and if so, return hit data.
        /// </summary>
        /// <param name="ray">Ray to check</param>
        /// <returns>Hit data (or null if no intersection)</returns>
        public RayHit Intersect(Ray ray)
        {
            // Implemented w/ Moller-Trumbore Algorithm
            // Adapted from https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection

            double epsilon = 0.000001;

            // Find the normal of the triangle
            Vector3 edge01 = this.v1 - this.v0;
            Vector3 edge02 = this.v2 - this.v0;

            Vector3 pvec = ray.Direction.Cross(edge02);
            double det = edge01.Dot(pvec);

            // Check if ray is parallel to the triangle
            if (Math.Abs(det) < epsilon) {
                return null;
            }

            // Compute u and v of the Barycentric coord
            double invDet = 1/det;

            Vector3 tvec = ray.Origin - this.v0;
            double u = tvec.Dot(pvec) * invDet;
            if (u < 0 || u > 1) return null;

            Vector3 qvec = tvec.Cross(edge01);
            double v = ray.Direction.Dot(qvec) * invDet;
            if (v < 0 || u + v > 1) return null;

            double t = edge02.Dot(qvec) * invDet;
            if (t <= 0) {
                //triangle is behind the Origin; not valid
                return null;    
            }

            Vector3 P = ray.Origin + t*ray.Direction;
            Vector3 normal = edge01.Cross(edge02).Normalized();
            
            return new RayHit(P, normal, ray.Direction, this.material);
        }
        

        /// <summary>
        /// Helper to check whether a point is in triangle
        /// using inside-outside test.
        /// </summary>
        /// <param name="P">Point to be checked</param>
        /// <param name="normal">Precomputed normal of the triangle</param>
        /// <returns>If point is in triangle</returns>
        private bool isInTriangle(Vector3 P, Vector3 normal) {
            Vector3 C; // perpendicular to triangle's plane

            // want C.Dot(normal) < 0

            // edge 0
            Vector3 edge01 = this.v1 - this.v0;
            Vector3 vP0 = P - this.v0;
            C = edge01.Cross(vP0);
            if (C.Dot(normal) < 0) {
                return false;
            }

            // edge 1
            Vector3 edge12 = this.v2-this.v1;
            Vector3 vP1 = P - this.v1;
            C = edge12.Cross(vP1);
            if (C.Dot(normal) < 0) {
                return false;
            }

            // edge 2
            Vector3 edge20 = this.v0-this.v2;
            Vector3 vP2 = P - this.v2;
            C = edge20.Cross(vP2);
            if (C.Dot(normal) < 0) {
                return false;
            }

            return true;

        }

        /// <summary>
        /// The material of the triangle.
        /// </summary>
        public Material Material { get { return this.material; } }
    }

        

}
