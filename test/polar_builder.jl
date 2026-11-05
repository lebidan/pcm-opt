# Run: julia --project=. test/polar_builder.jl
# py3gpp 0.6.0, nrCRCEncode(mask=0) + nrPolarEncode(E=128,nmax=7,iil=true).
# Row i is payload bit i set; hex packs codeword bits 0..127 MSB first.
using Test, MAT, LinearAlgebra
include(joinpath(@__DIR__, "..", "scripts", "build_polar_codes.jl"))

const BASIS_HEX = [
    "0312746503ed749afced8b9a03ed749a",
    "7d24b1e8421b8ed77d24b1e8421b8ed7",
    "fcb1ac1efc4eace1034eac1e03b1ace1",
    "42ee9fccbdee60ccbdee60cc42ee9fcc",
    "1ce68f7543b9d02a1ce68f7543b9d02a",
    "2d65783f2d65783f2d65783f2d65783f",
    "fc5fd172fca0d18d0caf21820c50217d",
    "5d57e316a2571c16d5df6b9e2adf949e",
    "509f11de5f901ed1509f11de5f901ed1",
    "9ccc78d79c3378285000b41b50ffb4e4",
    "90d759e1e7a02e9690d759e1e7a02e96",
    "875039dd875039dd875039dd875039dd",
    "3b2a2a3b3bd52ac491808091917f806e",
    "55f3336a66c0005955f3336a66c00059",
    "aa275a28aad85ad755d8a5d75527a528",
    "b5c162e9e09437bcb5c162e9e09437bc",
    "e5fe2c62e5fe2c62e5fe2c62e5fe2c62",
    "88e13a5388e13a5388e13a5388e13a53",
    "0c9642d8f39642270c9642d8f3964227",
    "fffa9369aafac669fffa9369aafac669",
]

function codeword_hex(bits)
    bytes = [sum(Int(bits[i+j]) << (7-j) for j in 0:7) for i in 1:8:length(bits)]
    join(string(b, base=16, pad=2) for b in bytes)
end

@testset "NR polar [128,20] CRC24C mother code" begin
    G, H, Ginv, perm = build()
    @test size(G) == (20, 128)
    @test size(H) == (108, 128)
    @test size(Ginv) == (128, 20)
    @test rank2(G) == 20
    @test rank2(H) == 108
    @test sum(H) == 2401
    @test !any(gfprod(G, transpose(H)))
    @test gfprod(G, Ginv) == Matrix{Bool}(I, 20, 20)
    @test [codeword_hex(G[i,:]) for i in 1:20] == BASIS_HEX
    # Include-based use must not rewrite the primary file; compare it to the build.
    saved = matread(joinpath(ROOT, "codes", "polar.128.20.crc24c.mat"))
    @test Int(saved["n"]) == 128 && Int(saved["k"]) == 20
    @test Matrix{Bool}(saved["G"]) == G
    @test Matrix{Bool}(saved["H"]) == H
    @test Matrix{Bool}(saved["Ginv"]) == Ginv
    @test vec(Int.(saved["input_permutation_0based"])) == perm .- 1
    check(G, H, Ginv, perm, "C")
end

const BASIS_HEX_40 = [
    "f49ecb51f49ecb51f49ecb51f49ecb51",
    "145be10e145b4ba4bef14ba4145b4ba4",
    "7149ee768eb6ee767149ee768eb6ee76",
    "217d994d217d994d217d994d217d994d",
    "ba1c4f49ef49b0b645e3b0b6ef49b0b6",
    "f82d98ed082d68edf82d98ed082d68ed",
    "742b39c6de8139c684db39c62e7139c6",
    "35d96925f9d9a52535d96925f9d9a525",
    "c11992ea6bb392ea0dd592eaa77f92ea",
    "9e50f19f34505b9f9e50f19f34505b9f",
    "0c32930da698930da698930d0c32930d",
    "dd343fd6779e3fd622cb3fd688613fd6",
    "c3015f3d3c01a03dc3015f3d3c01a03d",
    "0ab1657ea01b657efab1957e501b957e",
    "3806df41b8865fc13806df41b8865fc1",
    "a025597c0a8f597c6c25957cc68f957c",
    "a5989f0265585fc2a5989f0265585fc2",
    "fb647023fb647023fb647023fb647023",
    "3ef6076f945c076f94f6ad6f3e5cad6f",
    "a6c2e1850c68e18559c21e85f3681e85",
    "523c345af896345a92fcf49a3856f49a",
    "24dd2e77847d8ed724dd2e77847d8ed7",
    "be530c4114f90c411ef3ace1b459ace1",
    "36656f3c9ccf6f3cc6959fcc6c3f9fcc",
    "5a80a0daaa70502a5a80a0daaa70502a",
    "6588783f6588783f6588783f6588783f",
    "98e429f5324e29f5106ca17dbac6a17d",
    "53d9d852f973d8529f15149e35bf149e",
    "57b81659df309ed157b81659df309ed1",
    "c2b29e4e68189e4e681834e4c2b234e4",
    "afb7e25a637b2e96afb7e25a637b2e96",
    "cd56b9ddcd56b9ddcd56b9ddcd56b9dd",
    "6ca2ff91c608ff91935d006e39f7006e",
    "1f662af3b5cc80591f662af3b5cc8059",
    "98ca2528cd60252898ca2528cd602528",
    "94bfc8436b4037bc94bfc8436b4037bc",
    "bef02c62bef02c62bef02c62bef02c62",
    "d1b83a53d1b83a53d1b83a53d1b83a53",
    "5114c2270b4ec2275114c2270b4ec227",
    "ba954669dcf34669ba954669dcf34669",
]

const BASIS_HEX_72 = [
    "bb2e1e7bbbdeee8bbb2e1e7bbbdeee8b",
    "5e1b739c5114836c5e1b8c63aeeb7c93",
    "8d4b1478414bd878beb42787beb42787",
    "593b84e6f33b2ee60cc4d1190cc4d119",
    "e8835c371783a337e87c5cc8e87c5cc8",
    "a0b2825c9c4272aca0b2825c9c4272ac",
    "f0ebc37278634bfaf0ebc37278634bfa",
    "b4a2ee07b4a2ee07b4a2ee07b4a2ee07",
    "4366b396735683a67c598ca98ca97c59",
    "dd7878777dd8d8d78227272882272728",
    "9634444cccc4b4bc9634444cccc4b4bc",
    "cc6a222e00a6eee2cc6a222e00a6eee2",
    "8f8599938f8599938f8599938f859993",
    "15a24ff815a24ff81aad40f7ea5db007",
    "a0defa2e285672a6d7a98d59d7a98d59",
    "62649d9baea851575157aea85157aea8",
    "760c890c79fc79fc760c890c79fc79fc",
    "45a8efa8ef02450245a8efa8ef024502",
    "57b3678357b3678357b3678357b36783",
    "9bb594bac1efcee0cee0c1ef3e10311f",
    "69c2cccd963d333269c2cccd69c2cccd",
    "bff6b009bff6b009b0094ff6b0094ff6",
    "d83d18fde80de80dd83d18fde80de80d",
    "17def659e72e06a9242109a6d4d1f956",
    "4b19f50d1b4905fd4b19f50d1b4905fd",
    "f4c86d5104389da1a13792ae51c7625e",
    "dcd42381dcd4d371dcd42381dcd4d371",
    "059b509bf56ba06b0564af64f5945f94",
    "52d680ae2aae705e52d680ae2aae705e",
    "65f6b78e9a09487165f6b78e9a094871",
    "f5e8061b0518f6ebcad7f9e43a270914",
    "b1087d6e8d348d9eb1087d6e8d348d9e",
    "ee5ecb51ee5ecb51ee5ecb51ee5ecb51",
    "1bf444abeb04b45b44abbb54b45b4ba4",
    "1e2c1e864476ee761e2c1e864476ee76",
    "03bd994d03bd994d03bd994d03bd994d",
    "b0b6bfb940464f49bfb940464f49b0b6",
    "cde2981dc2ed68edcde2981dc2ed68ed",
    "633636c993c6c6391441c936e4b139c6",
    "43e9d5d53319a52543e9d5d53319a525",
    "ae7c9de55e8c6d159d4f621a6dbf92ea",
    "ce606b6ffe905b9fce606b6ffe905b9f",
    "63579c0293a76cf2360263fdc6f2930d",
    "e2a1c029e2a1c029e2a13fd6e2a13fd6",
    "ccf1f0cd9c01a03dccf1f0cd9c01a03d",
    "55d45a71a524aa816a2b658e9adb957e",
    "72b65f3172465fc172b65f3172465fc1",
    "a3403a7353b0ca83fcbf658c0c4f957c",
    "d7682732af985fc2d7682732af985fc2",
    "fda47023fda47023fda47023fda47023",
    "6153526091a3a2906eac5d9f9e5cad6f",
    "ee57697aee57697a99a81e8599a81e85",
    "0ba9c7650ba9c7653856f49a3856f49a",
    "188db227247d8ed7188db227247d8ed7",
    "db9609ee2b66f91e8e695c117e99ace1",
    "06009f3306009f3306ff9fcc06ff9fcc",
    "50800ada0a70502a50800ada0a70502a",
    "5548783f5548783f5548783f5548783f",
    "95492ef265b9de02ea36518d1ac6a17d",
    "60402ba160402ba15f7f149e5f7f149e",
    "70c091217f309ed170c091217f309ed1",
    "a7dd9b4b572d6bbbf882c414087234e4",
    "d9cb5ee6a9bb2e96d9cb5ee6a9bb2e96",
    "d096b9ddd096b9ddd096b9ddd096b9dd",
    "6608ff9196f80f616907f09e99f7006e",
    "4f3cb0697f0c80594f3cb0697f0c8059",
    "8027a2af70d7525ff750d5d807a02528",
    "f1d067eca18037bcf1d067eca18037bc",
    "98b02c6298b02c6298b02c6298b02c62",
    "d1b83a53d1b83a53d1b83a53d1b83a53",
    "688d01e4987df1145bbe32d7ab4ec227",
    "e366133ce366133cb6334669b6334669",
]

for (payload, edges, fixtures) in ((40, 2481, BASIS_HEX_40), (72, 2425, BASIS_HEX_72))
    @testset "NR polar [128,$payload] CRC24C mother code" begin
        G, H, Ginv, perm = build(; payload=payload)
        info, frozen = positions(payload)
        @test size(G) == (payload, 128)
        @test size(H) == (128-payload, 128)
        @test size(Ginv) == (128, payload)
        @test length(info) == payload + 24
        @test rank2(G) == payload
        @test rank2(H) == 128-payload
        @test sum(H) == edges
        @test !any(gfprod(G, transpose(H)))
        @test gfprod(G, Ginv) == Matrix{Bool}(I, payload, payload)
        @test [codeword_hex(G[i,:]) for i in 1:payload] == fixtures
        saved = matread(joinpath(ROOT, "codes", "polar.128.$payload.crc24c.mat"))
        @test Int(saved["n"]) == 128 && Int(saved["k"]) == payload
        @test Matrix{Bool}(saved["G"]) == G
        @test Matrix{Bool}(saved["H"]) == H
        @test Matrix{Bool}(saved["Ginv"]) == Ginv
        @test vec(Int.(saved["input_permutation_0based"])) == perm .- 1
        @test vec(Int.(saved["info_positions_0based"])) == info .- 1
        @test vec(Int.(saved["frozen_positions_0based"])) == frozen .- 1
        check(G, H, Ginv, perm, "C")
    end
end
@test_throws ArgumentError build(; payload=41)
