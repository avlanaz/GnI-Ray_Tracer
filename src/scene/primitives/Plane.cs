using System;

namespace RayTracer
{
    /// <summary>
    /// Class to represent an (infinite) plane in a scene.
    /// </summary>
    public class Plane : SceneEntity
    {
        private Vector3 center;
        private Vector3 normal;
        private Material material;

        /// <summary>
        /// Construct an infinite plane object.
        /// </summary>
        /// <param name="center">Position of the center of the plane</param>
        /// <param name="normal">Direction that the plane faces</param>
        /// <param name="material">Material assigned to the plane</param>
        public Plane(Vector3 center, Vector3 normal, Material material)
        {
            this.center = center;
            this.normal = normal.Normalized();
            this.material = material;
        }

        /// <summary>
        /// Determine if a ray intersects with the plane, and if so, return hit data.
        /// </summary>
        /// <param name="ray">Ray to check</param>
        /// <returns>Hit data (or null if no intersection)</returns>
        public RayHit Intersect(Ray ray)
        {
            // Code Written V
            double epsilon = 0.000001;

            // Any point on ray is given by P = o + t*d
            this.normal = this.normal.Normalized();

            // Check if plane normal parallel to ray (dir.Dot(N) ~ 0)
            double directionDotNormal = ray.Direction.Dot(this.normal);
            if (Math.Abs(directionDotNormal) < epsilon) {
                return null;
            }
            // To find t, solve t = ((o-p).Dot(N))/(d.Dot(N))
            double t = (this.center - ray.Origin).Dot(this.normal) / directionDotNormal;

            // If t < 0, plane is behind the Origin; not valid
            if (t <= 0) {
                return null;    
            }

            return new RayHit(ray.Origin + t*ray.Direction, this.normal, ray.Direction, this.material);
        }

        /// <summary>
        /// The material of the plane.
        /// </summary>
        public Material Material { get { return this.material; } }
    }

}
