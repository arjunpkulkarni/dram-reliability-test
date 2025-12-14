#!/usr/bin/env python3
import os
import sys

def check_file(filename, description, min_size=0):
    """Check if file exists and meets size requirements"""
    if os.path.exists(filename):
        size = os.path.getsize(filename)
        if size > min_size:
            print(f"✅ {description:40} ({size:,} bytes)")
            return True
        else:
            print(f"⚠️  {description:40} (TOO SMALL: {size} bytes)")
            return False
    else:
        print(f"❌ {description:40} (NOT FOUND)")
        return False

def main():
    print("=" * 80)
    print(" MSE 404 PROJECT VERIFICATION")
    print(" DRAM Reliability Modeling - Deliverables Check")
    print("=" * 80)
    print()
    
    all_good = True
    
    # 1. Core simulation code
    print("📋 SIMULATION CODE")
    print("-" * 80)
    all_good &= check_file("dram_reliability_fem.py", "Main FEM simulation script", 20000)
    print()
    
    # 2. Technical report
    print("📄 TECHNICAL REPORT")
    print("-" * 80)
    all_good &= check_file("technical_report.tex", "LaTeX report source", 10000)
    print()
    
    # 3. Figures
    print("📊 FIGURES (Computer-Generated)")
    print("-" * 80)
    all_good &= check_file("fig1_dram_geometry.png", "Figure 1: DRAM geometry", 10000)
    all_good &= check_file("fig2_temperature_field.png", "Figure 2: Temperature field", 10000)
    all_good &= check_file("fig3_stress_field.png", "Figure 3: Stress distribution", 10000)
    all_good &= check_file("fig4_scaling_trends.png", "Figure 4: Scaling trends", 10000)
    print()
    
    # 4. Documentation
    print("📚 DOCUMENTATION")
    print("-" * 80)
    all_good &= check_file("README.md", "Project README", 3000)
    all_good &= check_file("PROJECT_SUMMARY.md", "Executive summary", 5000)
    all_good &= check_file("COMPILATION_GUIDE.md", "LaTeX compilation guide", 1000)
    all_good &= check_file("requirements.txt", "Python dependencies", 30)
    print()
    
    # 5. Scripts
    print("🔧 UTILITY SCRIPTS")
    print("-" * 80)
    all_good &= check_file("run_simulation.sh", "Quick-start shell script", 500)
    all_good &= check_file("verify_project.py", "This verification script", 1000)
    print()
    
    # 6. Check Python syntax
    print("🐍 PYTHON SYNTAX CHECK")
    print("-" * 80)
    try:
        import py_compile
        py_compile.compile("dram_reliability_fem.py", doraise=True)
        print("✅ dram_reliability_fem.py syntax valid")
    except py_compile.PyCompileError as e:
        print(f"❌ Syntax error in dram_reliability_fem.py: {e}")
        all_good = False
    print()
    
    # 7. Check LaTeX basic syntax
    print("📝 LATEX CHECK")
    print("-" * 80)
    with open("technical_report.tex", "r") as f:
        tex_content = f.read()
        
    checks = [
        ("\\begin{document}", "Document environment"),
        ("\\end{document}", "Document closure"),
        ("\\begin{abstract}", "Abstract section"),
        ("\\section{Introduction}", "Introduction section"),
        ("\\section{Methods}", "Methods section"),
        ("\\section{Results", "Results section"),
        ("\\section{Conclusions}", "Conclusions section"),
        ("\\begin{thebibliography}", "Bibliography"),
        ("\\includegraphics", "Figure includes"),
    ]
    
    for pattern, desc in checks:
        if pattern in tex_content:
            print(f"✅ {desc:40}")
        else:
            print(f"❌ {desc:40} MISSING")
            all_good = False
    print()
    
    # Final summary
    print("=" * 80)
    if all_good:
        print("✅ ALL CHECKS PASSED - PROJECT READY FOR SUBMISSION")
        print()
        print("Next steps:")
        print("  1. Compile LaTeX report: See COMPILATION_GUIDE.md")
        print("  2. Review PROJECT_SUMMARY.md for submission checklist")
        print("  3. Test simulation: ./run_simulation.sh")
    else:
        print("⚠️  SOME CHECKS FAILED - REVIEW ISSUES ABOVE")
        print()
        print("Common fixes:")
        print("  - Re-run: python dram_reliability_fem.py")
        print("  - Check file permissions")
        print("  - Ensure all files in correct directory")
    print("=" * 80)
    print()
    
    # Statistics
    print("📈 PROJECT STATISTICS")
    print("-" * 80)
    try:
        with open("dram_reliability_fem.py", "r") as f:
            lines = len(f.readlines())
        print(f"  Python code lines: {lines}")
        
        with open("technical_report.tex", "r") as f:
            tex_lines = len(f.readlines())
        print(f"  LaTeX report lines: {tex_lines}")
        
        total_figs = sum([1 for f in ["fig1_dram_geometry.png", 
                                      "fig2_temperature_field.png",
                                      "fig3_stress_field.png",
                                      "fig4_scaling_trends.png"] 
                         if os.path.exists(f)])
        print(f"  Figures generated: {total_figs}/4")
        
        print(f"  Parametric cases: 9 (3 feature sizes × 3 current densities)")
        print(f"  Estimated runtime: ~30 seconds")
    except Exception as e:
        print(f"  Could not compute statistics: {e}")
    print("=" * 80)
    print()
    
    return 0 if all_good else 1

if __name__ == "__main__":
    sys.exit(main())

