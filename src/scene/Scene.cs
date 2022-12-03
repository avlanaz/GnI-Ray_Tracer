using System;
using System.Collections.Generic;

namespace RayTracer
{
    /// <summary>
    /// Class to represent a ray traced scene, including the objects,
    /// light sources, and associated rendering logic.
    /// </summary>
    public class Scene
    {
        private const int RecursionLimit = 4; // Limit for firing recursive rays
        private const int NAmbientSamples = 3;
        private const double Bias = 1e-4;   // Prevent hit with the same object
        private SceneOptions options;
        private ISet<SceneEntity> entities;
        private ISet<PointLight> lights;

        /// <summary>
        /// Construct a new scene with provided options.
        /// </summary>
        /// <param name="options">Options data</param>
        public Scene(SceneOptions options = new SceneOptions())
        {
            this.options = options;
            this.entities = new HashSet<SceneEntity>();
            this.lights = new HashSet<PointLight>();
        }

        /// <summary>
        /// Add an entity to the scene that should be rendered.
        /// </summary>
        /// <param name="entity">Entity object</param>
        public void AddEntity(SceneEntity entity)
        {
            this.entities.Add(entity);
        }

        /// <summary>
        /// Add a point light to the scene that should be computed.
        /// </summary>
        /// <param name="light">Light structure</param>
        public void AddPointLight(PointLight light)
        {
            this.lights.Add(light);
        }

        /// <summary>
        /// Render the scene to an output image. This is where the bulk
        /// of your ray tracing logic should go... though you may wish to
        /// break it down into multiple functions as it gets more complex!
        /// </summary>
        /// <param name="outputImage">Image to store render output</param>
        public void Render(Image outputImage)
        {
            // dotnet run -- -f tests/final_scene.txt -o output.png -h 100 -w 100 -x 2 -l
            DateTime start = DateTime.Now;
            double aspectRatio = outputImage.Width / outputImage.Height;
            int xLength = outputImage.Width*options.AAMultiplier;
            int yLength = outputImage.Height*options.AAMultiplier;
            Color[,] pixelColor = new Color[xLength, yLength];
            for (int y = 0; y < yLength; y++) {
                DateTime localStart = DateTime.Now;
                Console.WriteLine("y = " +y);
                for (int x = 0; x < xLength; x++) {
                    Ray pixelRay = ConstructPixelRay(x, y, aspectRatio, xLength, yLength);
                    RayHit hit = NearestHit(pixelRay);
                    if (hit == null) {
                        // Ray doesn't hit anything, move on to next ray
                        continue;
                    }
                    pixelColor[x,y] = ColorOfHit(hit, 0, x, y);
                }
                
                DateTime localEnd = DateTime.Now;
                TimeSpan localTs= localEnd-localStart;

                Console.WriteLine(localTs.TotalSeconds);  
            }

            // Averaging pixel colors for Anti-Aliasing
            for (int y = 0; y < outputImage.Height; y++) {
                for (int x = 0; x < outputImage.Width; x++) {
                    Color sumColor = new Color(0,0,0);

                    for (int j = 0; j < options.AAMultiplier; j++) {
                        for (int i = 0; i < options.AAMultiplier; i++) {
                            sumColor = sumColor + pixelColor[options.AAMultiplier*x+i, options.AAMultiplier*y+j];
                        }
                    }
                    Color avgColor = sumColor/(options.AAMultiplier*options.AAMultiplier);
                    outputImage.SetPixel(x, y, avgColor);
                }
            }
            DateTime end = DateTime.Now;
            TimeSpan ts= end-start;

            Console.WriteLine(ts.TotalSeconds);  
        }

        /// <summary>
        /// Helper to construct a ray for a certain pixel
        /// </summary>
        /// <param name="x">x location of pixel</param>
        /// <param name="y">y location of pixel</param>
        /// <param name="aspectRatio">precomputed aspect ratio of output image</param>
        /// <param name="outputImage">output image</param>
        /// <returns>Ray shoot through the pixel</returns>
        private Ray ConstructPixelRay(int x, int y, double aspectRatio, int width, int height) {
            //Construct a Ray for the current pixel
            double pixelX = (x+0.5)/width;
            double pixelY = (y+0.5)/height;

            double xPos = (pixelX*2) - 1;
            double yPos = 1 - (pixelY*2);
            double zPos = this.options.CameraAxis.Z;

            //Adjust to the 60deg FOV accordingly (using 60/2 for calculation)
            xPos = xPos*Math.Tan(30*Math.PI/180);
            yPos = yPos*Math.Tan(30*Math.PI/180) / aspectRatio;

            //Use camera position as Origin
            //Use the vector from Origin to Pixel Position as Direction
            Vector3 direction = (new Vector3(xPos,yPos,zPos) - this.options.CameraPosition);
            return new Ray(this.options.CameraPosition, direction);
        }

        /// <summary>
        /// Helper to determine the nearest object that the ray hits
        /// </summary>
        /// <param name="ray">the fired ray</param>
        /// <returns>The hit to the nearest object</returns>
        private RayHit NearestHit(Ray ray) {
            // Check hits for each object
            double nearestDist = Double.MaxValue;
            SceneEntity nearestEntity = null;
            RayHit hit = null;

            foreach(SceneEntity entity in this.entities) {
                RayHit thisHit = entity.Intersect(ray);
                if (thisHit != null) {
                    // It hits
                    double dist = (thisHit.Position-ray.Origin).LengthSq();
                    if (dist < nearestDist) {
                        nearestDist = dist;
                        hit = thisHit;
                        nearestEntity = entity;
                    }
                }
            }
            return hit;
        }

        /// <summary>
        /// Helper to handle color determination depending of the material
        /// </summary>
        /// <param name="hit">ray hit parameters</param>
        /// <param name="rayCounter">count how many reflections have been done</param>
        /// <returns>Color of the hit</returns>
        private Color ColorOfHit(RayHit hit, int rayCounter, int x, int y) {  
            // Diffuse handler
            if (hit.Material.Type == Material.MaterialType.Diffuse) {
                return DiffuseColor(hit, rayCounter, x, y);
            } 

            // Reflective handler
            else if (hit.Material.Type == Material.MaterialType.Reflective) {
                return ReflectiveColor(hit, rayCounter, x, y);
            } 

            // Refractive handler
            else if (hit.Material.Type == Material.MaterialType.Refractive) {
                return refractiveColor(hit, rayCounter, x, y);
            }

            // No recognised type, just use material color
            return hit.Material.Color;

        }
        
        /// <summary>
        /// Helper to determine color for diffuse material hit and ambient lighting
        /// Adapted from  https://www.scratchapixel.com/lessons/3d-basic-rendering/global-illumination-path-tracing/global-illumination-path-tracing-practical-implementation
        /// </summary>
        /// <param name="hit">ray hit parameters</param>
        /// <returns>Color of the hit</returns>
        private Color DiffuseColor(RayHit hit, int rayCounter, int x, int y) {
            if (rayCounter > RecursionLimit) {
                return new Color(0,0,0);
            }
            rayCounter+=1;
            Color directColor = DirectDiffuseColor(hit);
            Color indirectColor = new Color(0,0,0);

            // Return normal diffuse color if AL not enabled
            if (!(options.AmbientLightingEnabled)) {
                return directColor;
            }

            // Ambient Lighting Section
            var rand = new Random();
            Vector3[] coordVecs = ConstructNormalRelativeCoordinate(hit.Normal);

            if (rayCounter <= RecursionLimit) {
                // Fire sample rays for Ambient Lighting
                for (int i = 0; i < NAmbientSamples; i++) {
                    double e1 = rand.NextDouble();
                    double e2 = rand.NextDouble();

                    Vector3 sampleDir = UniformSampleHemisphereDirection(e1, e2);
                    Vector3 sampleAdjustedDir = new Vector3(
                        sampleDir.X * coordVecs[0].X + sampleDir.Y * coordVecs[1].X + sampleDir.Z * coordVecs[2].X, 
                        sampleDir.X * coordVecs[0].Y + sampleDir.Y * coordVecs[1].Y + sampleDir.Z * coordVecs[2].Y, 
                        sampleDir.X * coordVecs[0].Z + sampleDir.Y * coordVecs[1].Z + sampleDir.Z * coordVecs[2].Z);

                    Ray ambientRay = new Ray(hit.Position + sampleAdjustedDir*Bias, sampleAdjustedDir);
                    RayHit ambientRayHit = NearestHit(ambientRay);
                    Color ambientColor = new Color(0,0,0);
                    if (ambientRayHit != null) {
                        // multiply it by the PDF
                        ambientColor = ColorOfHit(ambientRayHit, rayCounter, x, y)*e1;
                    }
                    indirectColor += ambientColor;
                }
            }
            // Average the ambient colors
            indirectColor /= (float)NAmbientSamples;
            
            return (directColor/(Math.PI) + indirectColor*2)*(hit.Material.Color);
        }

        /// <summary>
        /// Helper to handle direct hit diffuse color
        /// </summary>
        /// <param name="hit">ray hit parameters</param>
        /// <param name="rayCounter">count how many reflections have been done</param>
        /// <returns>Color of the hit</returns>
        private Color DirectDiffuseColor(RayHit hit) {
            Color finalC = new Color(0.0f, 0.0f, 0.0f);

            Vector3 n = hit.Normal;
            if (n.Dot(hit.Incident) > 0) {
                // Hit from inside/backside, flip the normal
                n = -n;
            }

            foreach (PointLight light in this.lights) {
                // Check whether hit point is a shadow with respect to this light source
                bool hitSomething = false;
                Vector3 hitToLightPos = light.Position - hit.Position;
                double lightDist = hitToLightPos.LengthSq();

                // Fire a ray from hit point
                Ray hitRay = new Ray(hit.Position + n*Bias, hitToLightPos);

                // Check whether it hits an object nearer to a light source
                foreach(SceneEntity entity in this.entities) {
                    RayHit thisHit = entity.Intersect(hitRay);
                    if (thisHit != null) {
                        double entityDist = (hit.Position - thisHit.Position).LengthSq();
                        if (entityDist < lightDist) {
                            
                            // hits another entity closer to light
                            hitSomething = true;
                            break;
                        }
                    }
                }
                if (!hitSomething) {
                    double NDotL = n.Dot(hitToLightPos.Normalized());
                    Color C = light.Color * hit.Material.Color * Math.Max(0.0f, NDotL);
                    finalC = finalC + C;
                }
            }
            return finalC;
        }

        /// <summary>
        /// Helper to construct relative xyz coordinate with y+ being the normal
        /// </summary>
        /// <param name="N">ray normal</param>
        /// <returns>Basis vectors of the coordinate system</returns>
        private Vector3[] ConstructNormalRelativeCoordinate(Vector3 N) {
            // Create coordinate system relative to normal (DO THIS OUTSIDE S.T. DOESNT NEED TO RECOMPUTE)
            Vector3 Nt, Nb;
            if (Math.Abs(N.X) > Math.Abs(N.Y)) {
                Nt = (new Vector3(N.Z, 0.0f, -N.X))/Math.Sqrt(N.X*N.X + N.Z*N.Z);
            } else {
                Nt = (new Vector3(0, -N.Z, N.Y)) / Math.Sqrt(N.Y*N.Y + N.Z*N.Z);
            }
            Nb = N.Cross(Nt);

            return new Vector3[] {Nt, N, Nb};
        }

        /// <summary>
        /// Helper to construct relative xyz coordinate with y+ being the normal
        /// </summary>
        /// <param name="e1">random variable 1</param>
        /// <param name="e2">random variable 2</param>
        /// <returns>Random direction vector relative to a unit hemisphere</returns>
        private Vector3 UniformSampleHemisphereDirection(double e1, double e2) {
            // Construct random ray direction
            var rand = new Random();

            double sinTheta = Math.Sqrt(1 - e1 * e1); 
            double phi = 2 * Math.PI * e2; 
            double x = sinTheta * Math.Cos(phi); 
            double z = sinTheta * Math.Sin(phi); 
            return new Vector3(x, e1, z); // Hemisphere relative direction
        }

        /// <summary>
        /// Helper to determine color for reflective material hit
        /// Adapted from https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/reflection-refraction-fresnel
        /// </summary>
        /// <param name="hit">ray hit parameters</param>
        /// <param name="counter">count how many reflections have been done</param>
        /// <returns>Color of the hit</returns>
        private Color ReflectiveColor(RayHit hit, int rayCounter, int x, int y) {
            if (rayCounter > RecursionLimit) {
                return new Color(0.0f, 0.0f, 0.0f);
            }
            rayCounter += 1;
            
            // Reflection ray
            Vector3 reflectDir = hit.Incident - 2*hit.Normal*(hit.Normal.Dot(hit.Incident));
            Ray reflectRay = new Ray(hit.Position + Bias*reflectDir, reflectDir);

            // Get the nearest hits
            RayHit reflectHit = NearestHit(reflectRay);

            if (reflectHit == null) {
                // Doesnt hit anything else
                return new Color(0.0f, 0.0f, 0.0f);
            }

            // Material handling
            return ColorOfHit(reflectHit, rayCounter, x, y);
        }

        /// <summary>
        /// Helper to determine color for refractive material hit
        /// Adapted from https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/reflection-refraction-fresnel
        /// </summary>
        /// <param name="hit">ray hit parameters</param>
        /// <param name="counter">count how many reflections have been done</param>
        /// <returns>Color of the hit</returns>
        private Color refractiveColor(RayHit hit, int rayCounter, int x, int y) {
            if (rayCounter > RecursionLimit) {
                return new Color(0.0f, 0.0f, 0.0f);
            }
            rayCounter += 1;

            // Compute Fresnel Ratio
            double Fr = ComputeFresnel(hit);
            bool fromOutside = hit.Incident.Dot(hit.Normal) < 0;
            Vector3 bias = Bias*hit.Normal;

            // Compute Refraction Ray color of the Refraction Hit, if it isn't total internal reflection
            Color refractColor = new Color(0,0,0);
            if (Fr < 1) {
                double cosIn = Math.Clamp(hit.Incident.Dot(hit.Normal), -1, 1); 
                double etaIn = 1; 
                double etaOut = hit.Material.RefractiveIndex; 
                Vector3 n = hit.Normal; 
                if (cosIn < 0) { 
                    cosIn = -cosIn; // From outside
                }
                else { 
                    // From inside
                    etaIn = etaOut;
                    etaOut = 1;
                    n = -hit.Normal; 
                    } 
                double eta = etaIn / etaOut; 
                double k = 1 - eta * eta * (1 - cosIn * cosIn); 
                Vector3 dir = eta * hit.Incident + (eta * cosIn - Math.Sqrt(k)) * n;
                Vector3 orig = fromOutside ? hit.Position - bias : hit.Position + bias; 
                Ray refractRay = new Ray(orig, dir); 
                if (x == 275 && y == 256) {
                    Console.WriteLine(refractRay.Origin + " BBB " + refractRay.Direction);
                    }
                RayHit refractHit = NearestHit(refractRay);
                if (refractHit != null) {
                    refractColor = ColorOfHit(refractHit, rayCounter, x, y);
                }
            }

            // Compute color from Reflection Ray of the Refraction Hit
            Vector3 reflectDir = (hit.Incident - 2 * hit.Incident.Dot(hit.Normal) * hit.Normal).Normalized();
            Vector3 reflectOrig = fromOutside ? hit.Position + bias : hit.Position - bias; 
            Color reflectColor = new Color(0,0,0);

            Ray reflectRay = new Ray(reflectOrig, reflectDir);
            RayHit reflectHit = NearestHit(reflectRay);
            if (reflectHit != null) {
                reflectColor = ColorOfHit(reflectHit, rayCounter, x, y);
            }

            // Return the Fresnel color of the reflection and refraction
            return reflectColor*Fr + refractColor*(1-Fr);
        }

        /// <summary>
        /// Helper to compute fresnel ratio
        /// Adapted from https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/reflection-refraction-fresnel
        /// </summary>
        /// <param name="hit">ray hit parameters</param>
        /// <returns>The fresnel ratio</returns>
        private double ComputeFresnel(RayHit hit) {
            double cosIn = Math.Clamp(hit.Incident.Dot(hit.Normal), -1, 1);
            // Assume hit from outside
            double etaIn  = 1;
            double etaOut = hit.Material.RefractiveIndex;

            if (cosIn > 0) {
                // It's from inside
                etaIn = etaOut;
                etaOut = 1;
            }

            // Snell's law => sinOut
            double sinOut = etaIn/etaOut * Math.Sqrt(Math.Max(0f, 1-cosIn*cosIn));
            
            // Check if it's total internal reflection
            if (sinOut >= 1) {
                return 1;
            } 

            double cosOut = Math.Sqrt(Math.Max(0f, 1-sinOut*sinOut));
            cosIn = Math.Abs(cosIn);
            double Rs = ((etaOut * cosIn) - (etaIn * cosOut)) / ((etaOut * cosIn) + (etaIn * cosOut)); 
            double Rp = ((etaIn * cosIn) - (etaOut * cosOut)) / ((etaIn * cosIn) + (etaOut * cosOut));

            return (Rs*Rs + Rp*Rp)/2;
        }

    }
}
