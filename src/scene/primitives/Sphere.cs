using System;

namespace RayTracer
{
    /// <summary>
    /// Class to represent an (infinite) plane in a scene.
    /// </summary>
    public class Sphere : SceneEntity
    {
        private Vector3 center;
        private double radius;
        private Material material;

        /// <summary>
        /// Construct a sphere given its center point and a radius.
        /// </summary>
        /// <param name="center">Center of the sphere</param>
        /// <param name="radius">Radius of the spher</param>
        /// <param name="material">Material assigned to the sphere</param>
        public Sphere(Vector3 center, double radius, Material material)
        {
            this.center = center;
            this.radius = radius;
            this.material = material;
        }

        /// <summary>
        /// Determine if a ray intersects with the sphere, and if so, return hit data.
        /// </summary>
        /// <param name="ray">Ray to check</param>
        /// <returns>Hit data (or null if no intersection)</returns>
        public RayHit Intersect(Ray ray)
        {
            // Adapted from https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection
            // Using analytic solution

            Vector3 L = ray.Origin - this.center;
            double a = ray.Direction.Dot(ray.Direction);
            double b = 2 * ray.Direction.Dot(L);
            double c = L.Dot(L) - (this.radius * this.radius);

            // Solve the quadratic a*t^2 + b*t + c = 0
           
            double[] t = solveQuadratic(a,b,c);

            // Check if there's no intersection
            if (t == null) {
                return null;
            }
            

            // Use the t that gives the nearest intersection to origin
            double tHit;
            if (t[0] > t[1]) {
                double temp = t[0] ;
                t[0]  = t[1];
                t[1] = temp;
            }
 
            if (t[0]  < 0) { 
                 // t0 is negative, use t1
                t[0]  = t[1];
                if (t[0]  < 0) return null; // Both are negative, no valid intersection
            } 
 
            tHit = t[0];

            Vector3 point = ray.Origin + tHit*ray.Direction;
            Vector3 normal = (point-this.center).Normalized();

            return new RayHit(point, normal, ray.Direction, this.material);
        }

        /// <summary>
        /// Helper to solve the quadratic equation
        /// to find t for the intersection.
        /// </summary>
        /// <param name="a">component "a" of the quadratic</param>
        /// <param name="b">component "b" of the quadratic</param>
        /// <param name="b">component "c" of the quadratic</param>
        /// <returns>If point is in triangle</returns>
        private double[] solveQuadratic(double a, double b, double c) {
            double[] t = new double[2];
            double discriminant = b*b - 4*a*c;
            if (discriminant < 0) {
                // Ray doesn't intersect with sphere
                return null;
            } else if (discriminant == 0) {
                // Only one point of intersection
                t[0] = -0.5*b/a;
                t[1] = t[0];
            } else {
                // two points
                double q = (b > 0) ? -0.5 * (b + Math.Sqrt(discriminant)) : -0.5 * (b - Math.Sqrt(discriminant));
                t[0] = q/a;
                t[1] = c/q;
            }
            return t;
        }


        /// <summary>
        /// The material of the sphere.
        /// </summary>
        public Material Material { get { return this.material; } }
    }

}
