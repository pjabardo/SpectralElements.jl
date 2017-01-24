import Base.LinAlg: BlasInt, chkstride1
import Base.LinAlg.BLAS.@blasfunc
import Base.LinAlg.LAPACK.liblapack


# (GB) general banded matrices, LU decomposition and solver
for (pbtrf, pbtrs, elty) in
    ((:dpbtrf_,:dpbtrs_,:Float64),
     (:spbtrf_,:spbtrs_,:Float32),
     (:zpbtrf_,:zpbtrs_,:Complex128),
     (:cpbtrf_,:cpbtrs_,:Complex64))
    @eval begin
        # SUBROUTINE DPBTRF( UPLO, N, KD, AB, LDAB, INFO )
        # 
        # *       .. Scalar Arguments ..
        # *       CHARACTER          UPLO
        # *       INTEGER            INFO, KD, LDAB, N
        # *       ..
        # *       .. Array Arguments ..
        # *       DOUBLE PRECISION   AB( LDAB, * )
        function pbtrf!(uplo::Char, kd::Integer, AB::StridedMatrix{$elty})
            chkstride1(AB)
            n    = size(AB, 2)
            chkuplo(uplo)
            info = Ref{BlasInt}()
            ccall((@blasfunc($pbtrf), liblapack), Void,
                  (Ptr{UInt8}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt}, Ptr{BlasInt}),
                  &uplo, &n, &kd, AB, &max(1,stride(AB,2)), info)
            LAPACK.chklapackerror(info[])
            AB
        end

        # SUBROUTINE DPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
        # *
        # *       .. Scalar Arguments ..
        # *       CHARACTER          UPLO
        # *       INTEGER            INFO, KD, LDAB, LDB, N, NRHS
        # *       ..
        # *       .. Array Arguments ..
        # *       DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )
        function pbtrs!(trans::Char, kd::Integer, AB::StridedMatrix{$elty},
                        B::StridedVecOrMat{$elty})
            chkstride1(AB, B)
            LAPACK.chktrans(trans)
            info = Ref{BlasInt}()
            n    = size(AB,2)
            if n != size(B,1)
                throw(DimensionMismatch("Matrix AB has dimensions $(size(AB)), but right hand side matrix B has dimensions $(size(B))"))
            end
            ccall((@blasfunc($pbtrs), liblapack), Void,
                  (Ptr{UInt8}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, 
                   Ptr{$elty}, Ptr{BlasInt}, Ptr{$elty},   Ptr{BlasInt},
                   Ptr{BlasInt}),
                  &trans, &n, &kd, &size(B,2), AB, &max(1,stride(AB,2)), 
                  B, &max(1,stride(B,2)), info)
            LAPACK.chklapackerror(info[])
            B
        end
    end
end

"""
    pbtrf!(kl, ku, m, AB) -> (AB, ipiv)

Compute the LU factorization of a banded matrix `AB`. `kl` is the first
subdiagonal containing a nonzero band, `ku` is the last superdiagonal
containing one, and `m` is the first dimension of the matrix `AB`. Returns
the LU factorization in-place and `ipiv`, the vector of pivots used.
"""
pbtrf!(kl::Integer, ku::Integer, m::Integer, AB::StridedMatrix)

"""
    pbtrs!(trans, kl, ku, m, AB, ipiv, B)

Solve the equation `AB * X = B`. `trans` determines the orientation of `AB`. It may
be `N` (no transpose), `T` (transpose), or `C` (conjugate transpose). `kl` is the
first subdiagonal containing a nonzero band, `ku` is the last superdiagonal
containing one, and `m` is the first dimension of the matrix `AB`. `ipiv` is the vector
of pivots returned from `pbtrf!`. Returns the vector or matrix `X`, overwriting `B` in-place.
"""
pbtrs!(trans::Char, kl::Integer, ku::Integer, m::Integer, AB::StridedMatrix, ipiv::StridedVector{BlasInt}, B::StridedVecOrMat)

