## runtime benchmarks of megSAP


#### Note on DeepVariant GPU version
<!--- NA12878x2_93 --->
Calling variants with DeepVariant consists of three steps: `make_examples`, `call_variants`, `postprocess`.  
GPU acceleration only affects the `call_variants` step. Since `make_examples` takes up most of the overall runtime we don't support a DeepVariant GPU option for now.

Step-wise runtime comparison: GPU vs. CPU (WGS):

| Step           | GPU Time | CPU Time |
|----------------|------------|------------|
| `make_examples` | 377min19s | 466min 18s |
| `call_variants` | 10min 44s   | 43min 33s    |
| `postprocess`   | 2min 48s | 3min 29s        |
