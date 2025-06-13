"""
Test script for SUVTK commands
This script provides comprehensive testing of all suvtk commands with proper error handling
"""

import subprocess
# import sys
# import os
import argparse
from pathlib import Path
import time
from typing import List, Tuple, Optional

class Colors:
    """ANSI color codes for terminal output"""
    RED = '\033[0;31m'
    GREEN = '\033[0;32m'
    YELLOW = '\033[1;33m'
    BLUE = '\033[0;34m'
    MAGENTA = '\033[0;35m'
    CYAN = '\033[0;36m'
    NC = '\033[0m'  # No Color

class SUVTKTester:
    """Test runner for SUVTK commands"""
    
    def __init__(self, skip_commands: List[str] = None):
        self.input_dir = Path("test_examples/input")
        self.output_dir = Path("test_examples/output")
        self.database_dir = Path("test_examples/database")
        
        # Create directories
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.database_dir.mkdir(parents=True, exist_ok=True)
        
        # Test files
        self.test_fasta = self.input_dir / "test_sequences.fasta"
        self.single_fasta = self.input_dir / "YNP_partiti.fasta"
        
        # Results tracking
        self.results = []
        
        # Commands to skip
        self.skip_commands = skip_commands or []
        
    def run_command(self, cmd: List[str], description: str, 
                   capture_output: bool = True, 
                   timeout: int = 300) -> Tuple[bool, str, str]:
        """Run a command and return success status, stdout, stderr"""
        try:
            print(f"{Colors.BLUE}Running: {' '.join(cmd)}{Colors.NC}")
            
            result = subprocess.run(
                cmd, 
                capture_output=capture_output, 
                text=True, 
                timeout=timeout,
                cwd=Path.cwd()
            )
            
            success = result.returncode == 0
            return success, result.stdout, result.stderr
            
        except subprocess.TimeoutExpired:
            return False, "", f"Command timed out after {timeout} seconds"
        except Exception as e:
            return False, "", str(e)
    
    def test_command(self, cmd: List[str], description: str, 
                    expected_output_file: Optional[Path] = None) -> bool:
        """Test a single command and record results"""
        print(f"\n{Colors.YELLOW}Testing: {description}{Colors.NC}")
        print("-" * 60)
        
        start_time = time.time()
        success, stdout, stderr = self.run_command(cmd, description)
        end_time = time.time()
        
        duration = end_time - start_time
        
        if success:
            print(f"{Colors.GREEN}✓ SUCCESS{Colors.NC} ({duration:.2f}s)")
            if expected_output_file and expected_output_file.exists():
                size = expected_output_file.stat().st_size
                print(f"  Output file: {expected_output_file} ({size} bytes)")
        else:
            print(f"{Colors.RED}✗ FAILED{Colors.NC} ({duration:.2f}s)")
            if stderr:
                print(f"  Error: {stderr}")
        
        self.results.append({
            'description': description,
            'command': ' '.join(cmd),
            'success': success,
            'duration': duration,
            'stdout': stdout,
            'stderr': stderr
        })
        
        return success
    
    def test_help_commands(self):
        """Test all help commands to ensure they work"""
        if "help" in self.skip_commands:
            print(f"\n{Colors.YELLOW}=== Skipping Help Commands ==={Colors.NC}")
            return
            
        print(f"\n{Colors.MAGENTA}=== Testing Help Commands ==={Colors.NC}")
        
        help_commands = [
            (['suvtk', '--help'], "Main help"),
            (['suvtk', 'download-database', '--help'], "Download database help"),
            (['suvtk', 'taxonomy', '--help'], "Taxonomy help"),
            (['suvtk', 'features', '--help'], "Features help"),
            (['suvtk', 'virus-info', '--help'], "Virus info help"),
            (['suvtk', 'co-occurrence', '--help'], "Co-occurrence help"),
            (['suvtk', 'gbk2tbl', '--help'], "GBK2TBL help"),
            (['suvtk', 'comments', '--help'], "Comments help"),
            (['suvtk', 'table2asn', '--help'], "Table2ASN help"),
        ]
        
        for cmd, desc in help_commands:
            self.test_command(cmd, desc)
    
    def test_database_download(self):
        """Test database download"""
        if "database" in self.skip_commands or "download" in self.skip_commands:
            print(f"\n{Colors.YELLOW}=== Skipping Database Download cause it's slow ==={Colors.NC}")
            return
            
        print(f"\n{Colors.MAGENTA}=== Testing Database Download ==={Colors.NC}")
        
        cmd = ['suvtk', 'download-database', '-o', str(self.database_dir)]
        # Use a longer timeout for database download
        start_time = time.time()
        success, stdout, stderr = self.run_command(cmd, "Database download", timeout=600)
        end_time = time.time()
        
        duration = end_time - start_time
        
        if success:
            print(f"{Colors.GREEN}✓ SUCCESS{Colors.NC} ({duration:.2f}s)")
        else:
            print(f"{Colors.RED}✗ FAILED{Colors.NC} ({duration:.2f}s)")
            if stderr:
                print(f"  Error: {stderr}")
        
        self.results.append({
            'description': "Database download",
            'command': ' '.join(cmd),
            'success': success,
            'duration': duration,
            'stdout': stdout,
            'stderr': stderr
        })
    
    def test_main_commands(self):
        """Test main processing commands"""
        if "main" in self.skip_commands:
            print(f"\n{Colors.YELLOW}=== Skipping Main Commands ==={Colors.NC}")
            return
            
        print(f"\n{Colors.MAGENTA}=== Testing Main Commands ==={Colors.NC}")
        
        # Check if input files exist
        if not self.test_fasta.exists():
            print(f"{Colors.RED}Error: Test file not found: {self.test_fasta}{Colors.NC}")
            return
        
        if not self.single_fasta.exists():
            print(f"{Colors.RED}Error: Single test file not found: {self.single_fasta}{Colors.NC}")
            return
        
        # Test commands with main test file
        commands = [
            (['suvtk', 'taxonomy', '-i', str(self.test_fasta), '-o', 
              str(self.output_dir / 'taxonomy_results'), '-d', str(self.database_dir / 'suvtk_db')], 
             "Taxonomy assignment", self.output_dir / 'taxonomy_results', "taxonomy"),
            
            (['suvtk', 'features', '-i', str(self.test_fasta), '-o', 
              str(self.output_dir / 'features_results'), '-d', str(self.database_dir / 'suvtk_db')], 
             "Feature extraction", self.output_dir / 'features_results', "features"),
            
            (['suvtk', 'virus-info', '--taxonomy', str(self.output_dir / 'taxonomy_results' / 'taxonomy.tsv'), '-o', 
              str(self.output_dir / 'virus_info_results'), '-d', str(self.database_dir / 'suvtk_db')], 
             "Virus information", self.output_dir / 'virus_info_results', "virus-info"),
        ]
        
        for cmd, desc, expected_file, command_name in commands:
            if command_name not in self.skip_commands:
                self.test_command(cmd, desc, expected_file)
            else:
                print(f"{Colors.YELLOW}Skipping {desc}{Colors.NC}")
        
        # Test with single sequence file
        single_commands = [
            (['suvtk', 'taxonomy', '-i', str(self.single_fasta), '-o', 
              str(self.output_dir / 'single_taxonomy'), '-d', str(self.database_dir / 'suvtk_db')], 
             "Single sequence taxonomy", self.output_dir / 'single_taxonomy', "taxonomy"),
            
            (['suvtk', 'features', '-i', str(self.single_fasta), '-o', 
              str(self.output_dir / 'single_features'), '-d', str(self.database_dir / 'suvtk_db')], 
             "Single sequence features", self.output_dir / 'single_features', "features"),
        ]
        
        for cmd, desc, expected_file, command_name in single_commands:
            if command_name not in self.skip_commands:
                self.test_command(cmd, desc, expected_file)
            else:
                print(f"{Colors.YELLOW}Skipping {desc}{Colors.NC}")
    
    def create_test_files(self):
        """Create necessary test files for commands that require additional inputs"""
        
        # Create dummy MIUVIG file
        miuvig_path = self.output_dir / 'test_miuvig.tsv'
        miuvig_content = """MIUVIG_parameter\tvalue
source_uvig\tviral fraction metagenome (virome)
assembly_software\tmetaSPAdes;v3.15.3;kmer set 21,33,55,77, default otherwise
vir_ident_software\tgenomaD;1.7.0;score-calibration, default otherwise
size_frac\t0-0.8 um
virus_enrich_appr\tfiltration + centrifugation + DNAse + RNAse
nucl_acid_ext\t10.1038/srep16532
wga_amp_appr\tmda based"""
        
        with open(miuvig_path, 'w') as f:
            f.write(miuvig_content)
        
        # Create dummy assembly file
        assembly_path = self.output_dir / 'test_assembly.tsv'
        assembly_content = """Assembly_parameter\tvalue
StructuredCommentPrefix\tAssembly-Data
Assembly Method\tmetaSPAdes v. 3.15.3
Sequencing Technology\tIllumina NovaSeq 6000"""
        
        with open(assembly_path, 'w') as f:
            f.write(assembly_content)
        
        # Create simple abundance table for co-occurrence
        abundance_path = self.output_dir / 'test_abundance.tsv'
        abundance_content = """contig\tsample1\tsample2\tsample3
seq1\t10\t5\t0
seq2\t0\t15\t8
seq3\t5\t0\t12"""
        
        with open(abundance_path, 'w') as f:
            f.write(abundance_content)
        
        return miuvig_path, assembly_path, abundance_path
    
    def test_advanced_commands(self):
        """Test commands that require multiple input files"""
        if "advanced" in self.skip_commands:
            print(f"\n{Colors.YELLOW}=== Skipping Advanced Commands ==={Colors.NC}")
            return
            
        print(f"\n{Colors.MAGENTA}=== Testing Advanced Commands ==={Colors.NC}")
        
        # Create test files
        miuvig_path, assembly_path, abundance_path = self.create_test_files()
        
        # Check if taxonomy and features results exist
        taxonomy_path = self.output_dir / 'taxonomy_results' / 'miuvig_taxonomy.tsv'
        features_path = self.output_dir / 'features_results' / 'miuvig_features.tsv'
        
        if taxonomy_path.exists() and features_path.exists():
            # Test comments command
            if "comments" not in self.skip_commands:
                cmd = ['suvtk', 'comments', 
                       '-t', str(taxonomy_path),
                       '-f', str(features_path),
                       '-m', str(miuvig_path),
                       '-a', str(assembly_path),
                       '-o', str(self.output_dir / 'comments_results')]
                self.test_command(cmd, "Comments generation", 
                                self.output_dir / 'comments_results.cmt')
        else:
            print(f"{Colors.YELLOW}Skipping comments test - prerequisite files not found{Colors.NC}")
    
    def test_co_occurrence(self):
        """Test co-occurrence analysis (may need special handling)"""
        if "co-occurrence" in self.skip_commands:
            print(f"\n{Colors.YELLOW}=== Skipping Co-occurrence Analysis ==={Colors.NC}")
            return
            
        print(f"\n{Colors.MAGENTA}=== Testing Co-occurrence Analysis ==={Colors.NC}")
        
        # Create test abundance file
        _, _, abundance_path = self.create_test_files()
        
        # Test with abundance table instead of FASTA
        cmd = ['suvtk', 'co-occurrence', '-i', str(abundance_path), '-o', 
               str(self.output_dir / 'co_occurrence_results')]
        success = self.test_command(cmd, "Co-occurrence analysis", 
                                  self.output_dir / 'co_occurrence_results')
        
        if not success:
            print(f"{Colors.YELLOW}Note: Co-occurrence may need additional parameters{Colors.NC}")
    
    def print_summary(self):
        """Print test summary"""
        print(f"\n{Colors.MAGENTA}=== Test Summary ==={Colors.NC}")
        print("=" * 60)
        
        total_tests = len(self.results)
        successful_tests = sum(1 for r in self.results if r['success'])
        failed_tests = total_tests - successful_tests
        
        print(f"Total tests: {total_tests}")
        print(f"{Colors.GREEN}Successful: {successful_tests}{Colors.NC}")
        print(f"{Colors.RED}Failed: {failed_tests}{Colors.NC}")
        
        if failed_tests > 0:
            print(f"\n{Colors.RED}Failed Tests:{Colors.NC}")
            for result in self.results:
                if not result['success']:
                    print(f"  - {result['description']}")
                    if result['stderr']:
                        print(f"    Error: {result['stderr']}")
        
        print(f"\nOutput files in: {self.output_dir}")
        if self.output_dir.exists():
            output_files = list(self.output_dir.glob("*"))
            if output_files:
                print("Generated files:")
                for f in output_files:
                    size = f.stat().st_size if f.is_file() else 0
                    print(f"  - {f.name} ({size} bytes)")
    
    def run_all_tests(self):
        """Run complete test suite"""
        print(f"{Colors.CYAN}Starting SUVTK comprehensive test suite{Colors.NC}")
        print("=" * 60)
        
        # Test help commands first
        self.test_help_commands()
        
        # Test database download
        self.test_database_download()
        
        # Test main commands
        self.test_main_commands()
        
        # Test advanced commands (comments, etc.)
        self.test_advanced_commands()
        
        # Test co-occurrence separately
        self.test_co_occurrence()
        
        # Print summary
        self.print_summary()


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Comprehensive test suite for SUVTK commands",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python test_suvtk.py                                    # Run all tests
  python test_suvtk.py --skip-commands database          # Skip database download
  python test_suvtk.py --skip-commands database main     # Skip database and main commands
  python test_suvtk.py --skip-commands help taxonomy     # Skip help and taxonomy tests

Available commands to skip:
  help           - Skip all help command tests
  database       - Skip database download (also: download)
  download       - Skip database download (alias for database)
  main           - Skip all main processing commands
  advanced       - Skip advanced commands (comments, etc.)
  taxonomy       - Skip taxonomy assignment tests
  features       - Skip feature extraction tests
  virus-info     - Skip virus information tests
  comments       - Skip comments generation tests
  co-occurrence  - Skip co-occurrence analysis tests
        """)
    
    parser.add_argument(
        '--skip-commands',
        nargs='*',
        default=[],
        help='Commands to skip during testing. See examples below.'
    )
    
    args = parser.parse_args()
    
    tester = SUVTKTester(skip_commands=args.skip_commands)
    tester.run_all_tests()


if __name__ == "__main__":
    main() 