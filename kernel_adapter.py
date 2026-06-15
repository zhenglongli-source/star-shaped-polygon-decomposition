"""
Kernel adapter for the star-shaped decomposition project.

This version connects the star-shaped decomposition decomposition pipeline to the decomposition-aware
Lee-Preparata-inspired support-scan kernel computation.
"""

from kernel_support_scan import compute_kernel_support_scan, kernel_area


def compute_kernel(poly):
    kernel, _ = compute_kernel_support_scan(poly, return_trace=False)
    return kernel

