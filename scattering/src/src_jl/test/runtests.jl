using ScatteringI1
using HDF5
using Test
include("manual_correlators.jl")

function test_against_manually_constructed_correlation_matrix(infile)
    fid = h5open(infile)
     
    @testset "Test setup of correlation matrix elements" begin
        for ens in keys(fid)
            @testset "$ens" begin
                p_external = read(fid,"$ens/p_external")
                if "p(0,0,1)" ∈ p_external
                    ans_v1 = correlatorsp001(fid,ens;p=1)
                    ans_v2 = correlators_xyz(fid,ens;p=[0,0,1])
                    @test all(ans_v1 .== ans_v2)
                end
                if "p(0,0,2)" ∈ p_external
                    ans_v1 = correlatorsp001(fid,ens;p=2)
                    ans_v2 = correlators_xyz(fid,ens;p=[0,0,2])
                    @test all(ans_v1 .== ans_v2)
                end
                if "p(0,1,1)" ∈ p_external
                    ans_v1 = correlatorsp011(fid,ens;p=1)
                    ans_v2 = correlators_xyz(fid,ens;p=[0,1,1])
                    @test all(isapprox.(ans_v1,ans_v2))
                end
                if "p(1,1,0)" ∈ p_external
                    ans_v1 = correlatorsp110(fid,ens;p=1)
                    ans_v2 = correlators_xyz(fid,ens;p=[1,1,0])
                    @test all(isapprox.(ans_v1,ans_v2))
                end
            end
        end
    end
end

testfile = "../data/isospin1_merged.hdf5"
test_against_manually_constructed_correlation_matrix(testfile)