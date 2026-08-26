"""Length-tier -> shape_method ROUTING rule.

This is decided in Block 1 (so downstream capacity planning/reporting can use
it) but the `shape_method` column itself stays null until Block 2 actually
runs the corresponding fitter -- see SCHEMA.md "SHAPE LAYER". Basis:
results/rescue_characterization_summary.md (sub-read-length introns: classic
ratio suffices, AUC 0.955; 101-1000bp "awkward middle": two-point donor/
acceptor signal, AUC 0.677; >1000bp: full multi-bin shape, AUC 0.87-0.95).
"""
ROUTING = {
    "sub_read": "junction_only",
    "awkward_mid": "two_point",
    "long": "full_shape",
}


def route_shape_method(length_tier: str) -> str:
    return ROUTING[length_tier]
