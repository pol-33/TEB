# Mapper-Speed Profile Summary

## Run

| Field | Value |
| --- | --- |
| Date | 2026-04-24 19:37:03 CEST |
| Host | as02r5b15 |
| Slurm job id | 39610734 |
| Index | `./genome.dense.idx` |
| Reads | `../data/reads_1M.fastq` |
| Profile reads | 5000 |
| k | 1 |
| Mode | `all` |
| SIMD cap | `native` |
| Disable AVX512 | `0` |

Command template:

```bash
"/home/nct/nct01032/TEB/mapper-speed/mapper" -I "./genome.dense.idx" -i "./profile-results-dense/reads_subset.fastq" -o <tool-specific-output.sam> -k "1"
```

## Perf Stat

| Derived metric | Value |
| --- | --- |
| IPC | - |
| Branch miss rate | - |
| Cache miss rate | - |

| Event | Value | Unit |
| --- | --- | --- |
| task-clock:u | 1911.80 | msec |
| cycles:u | 3173459961 | - |
| instructions:u | 4162822933 | - |
| branches:u | 683092174 | - |
| branch-misses:u | 29638170 | - |
| cache-references:u | 37222010 | - |
| cache-misses:u | 29927560 | - |
| minor-faults:u | 892807 | - |
| major-faults:u | 0 | - |

## Perf Record

Raw files:

- `./profile-results-dense/perf-record/perf.data`
- `./profile-results-dense/perf-record/perf.report.txt`

### Top Hotspots

```text
    39.73%    38.90%  mapper   mapper                [.] mapper_speed::IndexView::offset_at
    14.34%     0.00%  mapper   [unknown]             [.] 0xffffffffa4000b40
    13.78%     0.07%  mapper   mapper                [.] mapper_speed::MapperEngine::search_orientation
    12.31%     0.66%  mapper   mapper                [.] mapper_speed::MapperEngine::prefilter_candidate
    11.55%    10.35%  mapper   mapper                [.] mapper_speed::banded_score_only
    11.09%    10.79%  mapper   mapper                [.] mapper_speed::IndexView::extract_sequence
     9.92%     9.59%  mapper   mapper                [.] mapper_speed::MapperEngine::generate_candidates
     7.09%     7.09%  mapper   mapper                [.] std::__introsort_loop<__gnu_cxx::__normal_iterator<mapper_speed::SeedAnchor*, std::vector<mapper_speed::SeedAnchor, std::allocator<mapper_speed::SeedAnchor> > >, long, __gnu_cxx::__ops::_Iter_comp_iter<mapper_speed::MapperEngine::generate_candidates(std::vector<mapper_speed::SeedSpec, std::allocator<mapper_speed::SeedSpec> > const&, unsigned long, int, std::vector<mapper_speed::CandidateInfo, std::allocator<mapper_speed::CandidateInfo> >&)::{lambda(mapper_speed::SeedAnchor const&, mapper_speed::SeedAnchor const&)#3}> >
     5.25%     5.14%  mapper   mapper                [.] mapper_speed::IndexView::occurrence_count
     2.95%     2.95%  mapper   mapper                [.] mapper_speed::MapperEngine::collect_seeds
     2.73%     2.73%  mapper   mapper                [.] std::__introsort_loop<__gnu_cxx::__normal_iterator<mapper_speed::MapperEngine::generate_candidates(std::vector<mapper_speed::SeedSpec, std::allocator<mapper_speed::SeedSpec> > const&, unsigned long, int, std::vector<mapper_speed::CandidateInfo, std::allocator<mapper_speed::CandidateInfo> >&)::ChainCluster*, std::vector<mapper_speed::MapperEngine::generate_candidates(std::vector<mapper_speed::SeedSpec, std::allocator<mapper_speed::SeedSpec> > const&, unsigned long, int, std::vector<mapper_speed::CandidateInfo, std::allocator<mapper_speed::CandidateInfo> >&)::ChainCluster, std::allocator<mapper_speed::MapperEngine::generate_candidates(std::vector<mapper_speed::SeedSpec, std::allocator<mapper_speed::SeedSpec> > const&, unsigned long, int, std::vector<mapper_speed::CandidateInfo, std::allocator<mapper_speed::CandidateInfo> >&)::ChainCluster> > >, long, __gnu_cxx::__ops::_Iter_comp_iter<mapper_speed::MapperEngine::generate_candidates(std::vector<mapper_speed::SeedSpec, std::allocator<mapper_speed::SeedSpec> > const&, unsigned long, int, std::vector<mapper_speed::CandidateInfo, std::allocator<mapper_speed::CandidateInfo> >&)::{lambda(mapper_speed::MapperEngine::generate_candidates(std::vector<mapper_speed::SeedSpec, std::allocator<mapper_speed::SeedSpec> > const&, unsigned long, int, std::vector<mapper_speed::CandidateInfo, std::allocator<mapper_speed::CandidateInfo> >&)::ChainCluster const&, mapper_speed::MapperEngine::generate_candidates(std::vector<mapper_speed::SeedSpec, std::allocator<mapper_speed::SeedSpec> > const&, unsigned long, int, std::vector<mapper_speed::CandidateInfo, std::allocator<mapper_speed::CandidateInfo> >&)::ChainCluster const&)#5}> >
     2.50%     2.50%  mapper   mapper                [.] mapper_speed::(anonymous namespace)::myers_popcnt
     1.35%     1.18%  mapper   mapper                [.] mapper_speed::IndexView::count_mismatches_packed
     1.31%     1.20%  mapper   mapper                [.] std::vector<int, std::allocator<int> >::_M_fill_assign
     1.27%     0.00%  mapper   [unknown]             [.] 0x00000000014341b0
     1.22%     1.22%  mapper   mapper                [.] std::__introsort_loop<__gnu_cxx::__normal_iterator<mapper_speed::SeedSpec*, std::vector<mapper_speed::SeedSpec, std::allocator<mapper_speed::SeedSpec> > >, long, __gnu_cxx::__ops::_Iter_comp_iter<mapper_speed::MapperEngine::collect_seeds(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, int, std::vector<mapper_speed::SeedSpec, std::allocator<mapper_speed::SeedSpec> >&)::{lambda(mapper_speed::SeedSpec const&, mapper_speed::SeedSpec const&)#1}> >
     1.02%     0.82%  mapper   mapper                [.] mapper_speed::MapperEngine::try_upfront_exactish
     0.93%     0.93%  mapper   [unknown]             [k] 0xffffffffa400108d
     0.58%     0.58%  mapper   mapper                [.] mapper_speed::IndexView::positions_for
     0.48%     0.48%  mapper   mapper                [.] mapper_speed::IndexView::chromosome_for_position
     0.44%     0.44%  mapper   libc.so.6             [.] __memchr_evex
     0.37%     0.27%  mapper   libstdc++.so.6.0.29   [.] std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >::_M_replace_aux
     0.25%     0.11%  mapper   mapper                [.] mapper_speed::IndexView::chromosome
     0.25%     0.25%  mapper   mapper                [.] std::__insertion_sort<__gnu_cxx::__normal_iterator<mapper_speed::SeedSpec*, std::vector<mapper_speed::SeedSpec, std::allocator<mapper_speed::SeedSpec> > >, __gnu_cxx::__ops::_Iter_comp_iter<mapper_speed::MapperEngine::collect_seeds(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, int, std::vector<mapper_speed::SeedSpec, std::allocator<mapper_speed::SeedSpec> >&)::{lambda(mapper_speed::SeedSpec const&, mapper_speed::SeedSpec const&)#2}> >
     0.23%     0.00%  mapper   [unknown]             [.] 0x0000002c00010003
```

