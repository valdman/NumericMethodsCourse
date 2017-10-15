using System;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using MathNet.Numerics.Properties;

namespace Lab1.QrHelpers
{
    public abstract class QR : QR<double>
    {
        protected QR(Matrix<double> q, Matrix<double> rFull, QRMethod method)
            : base(q, rFull, method)
        {
        }

        public override double Determinant
        {
            get
            {
                if (this.FullR.RowCount != this.FullR.ColumnCount)
                    throw new ArgumentException(Resources.ArgumentMatrixSquare);
                var num = 1.0;
                for (var index = 0; index < this.FullR.ColumnCount; ++index)
                {
                    num *= this.FullR.At(index, index);
                    if (Math.Abs(this.FullR.At(index, index)).AlmostEqual(0.0))
                        return 0.0;
                }
                return Math.Abs(num);
            }
        }

        public override bool IsFullRank
        {
            get
            {
                for (var index = 0; index < this.FullR.ColumnCount; ++index)
                {
                    if (Math.Abs(this.FullR.At(index, index)).AlmostEqual(0.0))
                        return false;
                }
                return true;
            }
        }
    }
}