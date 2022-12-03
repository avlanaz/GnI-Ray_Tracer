using System;

namespace RayTracer
{
    /// <summary>
    /// Immutable structure to represent a three-dimensional vector.
    /// </summary>
    public readonly struct Vector3
    {
        private readonly double x, y, z;

        /// <summary>
        /// Construct a three-dimensional vector.
        /// </summary>
        /// <param name="x">X component</param>
        /// <param name="y">Y component</param>
        /// <param name="z">Z component</param>
        public Vector3(double x, double y, double z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }

        /// <summary>
        /// Convert vector to a readable string.
        /// </summary>
        /// <returns>Vector as string in form (x, y, z)</returns>
        public override string ToString()
        {
            return "(" + this.x + "," + this.y + "," + this.z + ")";
        }

        /// <summary>
        /// Compute the length of the vector squared.
        /// This should be used if there is a way to perform a vector
        /// computation without needing the actual length, since
        /// a square root operation is expensive.
        /// </summary>
        /// <returns>Length of the vector squared</returns>
        public double LengthSq()
        {
            // Code Written V
            return (this.x*this.x) + (this.y*this.y) + (this.z*this.z);
        }

        /// <summary>
        /// Compute the length of the vector.
        /// </summary>
        /// <returns>Length of the vector</returns>
        public double Length()
        {
            // Code Written V
            return Math.Sqrt(this.LengthSq());
        }

        /// <summary>
        /// Compute a length 1 vector in the same direction.
        /// </summary>
        /// <returns>Normalized vector</returns>
        public Vector3 Normalized()
        {
            // Code Written V
            double xN = x/this.Length();
            double yN = y/this.Length();
            double zN = z/this.Length();

            // Can use Fast Inverse Square Root for faster algorithm
            return new Vector3(xN, yN, zN);
        }

        /// <summary>
        /// Compute the dot product with another vector.
        /// </summary>
        /// <param name="with">Vector to dot product with</param>
        /// <returns>Dot product result</returns>
        public double Dot(Vector3 with)
        {
            // Code Written V
            return (this.x*with.x) + (this.y*with.y) + (this.z*with.z);
        }

        /// <summary>
        /// Compute the cross product with another vector.
        /// </summary>
        /// <param name="with">Vector to cross product with</param>
        /// <returns>Cross product result</returns>
        public Vector3 Cross(Vector3 with)
        {
            // Code Written V
            double xCross = (this.y*with.z) - (this.z*with.y);
            double yCross = (this.x*with.z) - (this.z*with.x);
            double zCross = (this.x*with.y) - (this.y*with.x);
            return new Vector3(xCross, -yCross, zCross);
        }

        /// <summary>
        /// Sum two vectors together (using + operator).
        /// </summary>
        /// <param name="a">First vector</param>
        /// <param name="b">Second vector</param>
        /// <returns>Summed vector</returns>
        public static Vector3 operator +(Vector3 a, Vector3 b)
        {
            // Code Written V
            double xAdd = a.x + b.x;
            double yAdd = a.y + b.y;
            double zAdd = a.z + b.z;
            return new Vector3(xAdd, yAdd, zAdd);
        }

        /// <summary>
        /// Negate a vector (using - operator)
        /// </summary>
        /// <param name="a">Vector to negate</param>
        /// <returns>Negated vector</returns>
        public static Vector3 operator -(Vector3 a)
        {
            // Code Written V
            double xNeg = -a.x;
            double yNeg = -a.y;
            double zNeg = -a.z;
            return new Vector3(xNeg, yNeg, zNeg);
        }

        /// <summary>
        /// Subtract one vector from another.
        /// </summary>
        /// <param name="a">Original vector</param>
        /// <param name="b">Vector to subtract</param>
        /// <returns>Subtracted vector</returns>
        public static Vector3 operator -(Vector3 a, Vector3 b)
        {
            // Code Written V
            double xSub = a.x - b.x;
            double ySub = a.y - b.y;
            double zSub = a.z - b.z;
            return new Vector3(xSub, ySub, zSub);
        }

        /// <summary>
        /// Multiply a vector by a scalar value.
        /// </summary>
        /// <param name="a">Original vector</param>
        /// <param name="b">Scalar multiplier</param>
        /// <returns>Multiplied vector</returns>
        public static Vector3 operator *(Vector3 a, double b)
        {
            // Code Written V
            double xMult = b*a.x;
            double yMult = b*a.y;
            double zMult = b*a.z;
            return new Vector3(xMult, yMult, zMult);
        }

        /// <summary>
        /// Multiply a vector by a scalar value (opposite operands).
        /// </summary>
        /// <param name="b">Scalar multiplier</param>
        /// <param name="a">Original vector</param>
        /// <returns>Multiplied vector</returns>
        public static Vector3 operator *(double b, Vector3 a)
        {
            // Code Written V
            double xMult = b*a.x;
            double yMult = b*a.y;
            double zMult = b*a.z;
            return new Vector3(xMult, yMult, zMult);
        }

        /// <summary>
        /// Divide a vector by a scalar value.
        /// </summary>
        /// <param name="a">Original vector</param>
        /// <param name="b">Scalar divisor</param>
        /// <returns>Divided vector</returns>
        public static Vector3 operator /(Vector3 a, double b)
        {
            // Code Written V
            double xDiv = a.x/b;
            double yDiv = a.y/b;
            double zDiv = a.z/b;
            return new Vector3(xDiv, yDiv, zDiv);
        }

        /// <summary>
        /// X component of the vector.
        /// </summary>
        public double X { get { return this.x; } }

        /// <summary>
        /// Y component of the vector.
        /// </summary>
        public double Y { get { return this.y; } }

        /// <summary>
        /// Z component of the vector.
        /// </summary>
        public double Z { get { return this.z; } }
    }
}
