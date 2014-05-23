"""
OptMAGE

Originally released 5/4/2009 by:
    Harris Wang (hhwang@genetics.med.harvard.edu)
    Church Lab
    Harvard Medical School
    Boston, MA 02115, USA
    Copyright (C) 2009

Changes made in 2011-2012 by:
    Marc J. Lajoie (mlajoie@genetics.med.harvard.edu)
    Daniel B. Goodman (dbg@mit.edu)

Ported to python in 2013 by:
    Gleb Kuznetsov (gleb@mit.edu)
"""

import csv
import math
import os
import subprocess

import argparse

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq

PWD = os.path.dirname(os.path.realpath(__file__))
DATA = os.path.join(PWD, 'data')

# I/O related constants/defaults.
DEFAULT_GENOME = os.path.join(DATA, 'mg1655.fasta')
DEFAULT_OLIGO_TARGET_FILENAME = os.path.join(DATA, 'input_targets.csv')
DEFAULT_OLIGO_RESULT_OUTPUT_FILENAME = os.path.join(DATA, 'out_oligos.csv')

# Other defaults/constants.
DEFAULT_OLIGO_SIZE = 90
DEFAULT_MIN_SS_DG = -12.0
DEFAULT_MUT_LOC_MAX = 20
DEFAULT_NUM_PHOSPHOROTHIOATE = 2

# Replichore-related constants.
DEFAULT_REPLICATION_ORIGIN = (3923767, 3923998)
DEFAULT_REPLICATION_TERMINUS = (1588774, 1588801)

# NOTE: These were the coordinates used in the old script,
# but they appear to be wrong.
# DEFAULT_REPLICATION_ORIGIN = (3932974, 3933205)
# DEFAULT_REPLICATION_TERMINUS = (1597981, 1598008)

# Symbols
PHOSPHOROTHIOATE_SYMBOl = '*'


OPT_MAGE_MUTATION_TYPE__REFERENCE = 'R'
OPT_MAGE_MUTATION_TYPE__SUBSTITUTION = 'M'
OPT_MAGE_MUTATION_TYPE__DELETION = 'D'
OPT_MAGE_MUTATION_TYPE__INSERTION = 'I'

VALID_MUTATION_TYPES = set([
    OPT_MAGE_MUTATION_TYPE__REFERENCE,
    OPT_MAGE_MUTATION_TYPE__SUBSTITUTION,
    OPT_MAGE_MUTATION_TYPE__DELETION,
    OPT_MAGE_MUTATION_TYPE__INSERTION
])

# Map from acceptable values for the strand parameter to the representation
# the logic of the script uses.
STRAND_INPUT_INTERPRETATION_MAP = {
        '+': 1,
        '1': 1,
        '+1': 1,
        '-': -1,
        '-1': -1,
        1: 1,
        -1: -1
}


OLIGO_TARGET_REQUIRED_PARAMS = set([
    'target_id',
    'strand',
    'start',
    'end',
    'mutation_type',
])


class OptMAGEConfig(object):
    """Object that stores configuration for a design run.
    """

    def __init__(self,
            oligo_size=DEFAULT_OLIGO_SIZE,
            min_ss_dG=DEFAULT_MIN_SS_DG,
            oligo_end_buffer_distance=DEFAULT_MUT_LOC_MAX,
            num_phosphorothioate_bonds=DEFAULT_NUM_PHOSPHOROTHIOATE,
            auto_calc_replichore=True,
            ref_genome_source_location=DEFAULT_GENOME,
            replication_origin=DEFAULT_REPLICATION_ORIGIN,
            replication_terminus=DEFAULT_REPLICATION_TERMINUS):
        """Constructor."""
        # The size of the oligo. This is typically 90.
        self.oligo_size = oligo_size

        # Minimum free energy that we want to allow an oligo to have in order
        # to avoid strong secondary structure.
        self.min_ss_dG = min_ss_dG

        # The minimum number of base pairsa mutation and the end of
        # an oligo. Previous experiments have shown that there is a sharp drop
        # in replacement efficiency if the mutation gets within ~15 base pairs
        # of the oligo.
        self.oligo_end_buffer_distance = oligo_end_buffer_distance

        # The number of terminal 5' phosphorothioate bonds.
        self.num_phosphorothioate_bonds = num_phosphorothioate_bonds

        # 0 = False, 1 = True
        # If True, use the mutation position to automatically calculate
        # replichore information. Otherwise, use the replichore value specified
        # in the oligo target file.
        self.should_calc_replichore = auto_calc_replichore

        # The genome.
        self.set_genome_record_from_source(ref_genome_source_location)

        # Replichore-related.
        self.replication_origin = replication_origin
        self.replication_terminus = replication_terminus

    @classmethod
    def build_from_args(cls, args):
        """Factory method for creating a config object from parsed args.
        """
        return cls(
            oligo_size=args.oligo_size,
            min_ss_dG=args.min_ss_dG,
            oligo_end_buffer_distance=args.mut_loc_max,
            num_phosphorothioate_bonds=args.num_thio,
            auto_calc_replichore=(not args.manually_calc_replichore),
            ref_genome_source_location=args.ref_genome,
            replication_origin=args.replication_origin,
            replication_terminus=args.replication_terminus)

    def set_genome_record_from_source(self, genome_source):
        """Sets the genome record from source file path.
        """
        with open(genome_source) as genome_fh:
            self.genome_record = SeqIO.read(genome_fh, 'fasta')


class OligoTarget(object):
    """Object that specifies the properties of the Oligo to create.
    """

    MIN_EXPECTED_ARGS = 6

    def __init__(self, optMAGE_config, params):
        """Constructor.

        Args:
            optMAGE_config: An OptMAGEConfig object.
            params: Dictionary with keys:
                * target_id
                * strand
                * start # 1-indexed
                * end # 1-indexed (base after last base in seq)
                * mutation_type
                * replichore (optional)
                * mutation_seq (optional)

        Returns:
            An OligoTarget object instance.
        """
        assert not OLIGO_TARGET_REQUIRED_PARAMS - set(params.keys())
        self.target_id = params['target_id']

        self.strand = STRAND_INPUT_INTERPRETATION_MAP[params['strand']]
        assert self.strand == 1 or self.strand == -1

        # pythonic
        self.start = int(params['start']) - 1
        assert self.start >= 0

        # pythonic
        self.end = int(params['end']) - 1
        assert self.end <= len(optMAGE_config.genome_record)
        assert self.end >= self.start, "Bad input. End is less than start."

        self.mutation_type = params['mutation_type']
        assert self.mutation_type in VALID_MUTATION_TYPES, (
                "Invalid mutation type %s" % self.mutation_type)

        self.original_seq = optMAGE_config.genome_record.seq[
                self.start:self.end]

        # Determine the mutation sequence. If mutation_type is R, "from
        # reference genome", then we just take that slice from directly
        # from the reference genome. Otherwise we get it from the provided
        # argument.
        if self.mutation_type == OPT_MAGE_MUTATION_TYPE__REFERENCE:
            self.mutation_seq = optMAGE_config.genome_record.seq[
                    self.start:self.end]
            if self.strand == -1:
                self.mutation_seq = self.mutation_seq.reverse_complement()
        elif self.mutation_type == OPT_MAGE_MUTATION_TYPE__DELETION:
            self.mutation_seq = ''
        else:
            self.mutation_seq = Seq(params['mutation_seq'], generic_dna)
        # Basic validation.
        assert len(self.mutation_seq) / 2 <= \
                optMAGE_config.oligo_size / 2 - \
                        optMAGE_config.oligo_end_buffer_distance, \
                "Mutation is too large."

        # Set replichore. Must be either 1 or 2.
        if optMAGE_config.should_calc_replichore:
            self.replichore = self._calc_replichore(
                    optMAGE_config, self.start, self.end)
        else:
            self.replichore = int(params['replichore'])
        assert self.replichore == 1 or self.replichore == 2


    def _calc_replichore(self, optMAGE_config, start, end):
        """Calculate the replichore from start and end positions.
        """
        mutation_midpoint = (end - start) / 2 + start
        if (mutation_midpoint < optMAGE_config.replication_terminus[0] or
                mutation_midpoint > optMAGE_config.replication_origin[1]):
            return 1
        elif (mutation_midpoint > optMAGE_config.replication_terminus[1] and
                mutation_midpoint < optMAGE_config.replication_origin[0]):
            return 2
        else:
            # TODO: Look into literature and make an actual call here, e.g.:
            #     * http://www.ncbi.nlm.nih.gov/pubmed/2846183
            #     * http://www.annualreviews.org/doi/abs/10.1146/annurev.ge.26.120192.002311
            raise RuntimeError(
                    "Unable to determine replichore for (start, end) %s." %
                            ((start, end),))


    @classmethod
    def parse_input_and_create_list(cls, optMAGE_config, input_path):
        """Factory method that returns a list of OligoTarget objects.
        """
        oligo_target_list = []
        with open(input_path) as fh:
            fh.readline() # Skip the header.
            current_line = fh.readline()
            while current_line:
                # Grag the arg_list. Currently we support both space-separated
                # and comma-separated inputs.
                arg_list = current_line.split()
                if not len(arg_list) >= cls.MIN_EXPECTED_ARGS:
                    arg_list = current_line.split(',')
                    arg_list = [arg.strip() for arg in arg_list]
                    assert len(arg_list) >= len(OLIGO_TARGET_REQUIRED_PARAMS), (
                        "Invalid oligo target input.")

                # Create the params object from the arg_list.
                params = {
                    'target_id': arg_list[0],
                    'strand': arg_list[1],
                    'replichore': arg_list[2],
                    'start': arg_list[3],
                    'end': arg_list[4],
                    'mutation_type': arg_list[5]
                }
                if len(arg_list) > 6:
                    params['mutation_seq'] = arg_list[6]

                # Create the OligoTarget object.
                oligo_target_list.append(cls(optMAGE_config, params))

                current_line = fh.readline()

        # Return the complete list.
        return oligo_target_list


class OligoResult(object):
    """Object representing an oligo that has been designed.
    """

    def __init__(self, target_id, start, end, strand, replichore,
            mutation_type, ss_dG, oligo_size, oligo_seq, original_seq,
            mutation_seq, predicted_replacement_efficiency):
        """Constructor."""
        self.target_id = target_id
        self.start = start
        self.end = end
        self.strand = strand
        self.replichore = replichore
        self.mutation_type = mutation_type
        self.ss_dG = ss_dG
        self.oligo_size = oligo_size
        self.oligo_seq = oligo_seq
        self.original_seq = original_seq
        self.mutation_seq = mutation_seq
        self.predicted_replacement_efficiency = predicted_replacement_efficiency


class OligoGenerator(object):
    """Object that encapsulates the flow for generating oligos based
    on configuration and specifications.
    """

    def __init__(self, optMAGE_config):
        self.config = optMAGE_config


    def generate_oligo(self, oligo_target):
        """Generates an oligo given the config and target specs.

        Returns:
            An OligoResult object.
        """
        # Caculate the candidate block from which we'll select the oligo.
        block_seq = self.get_candidate_block_seq(oligo_target)

        # Determine the oligo cut-off from the block and report its final
        # free energy.
        oligo_seq_result, _ = self.determine_oligo_from_block(block_seq)
        oligo_seq = oligo_seq_result['oligo_seq']
        ss_dG = oligo_seq_result['ss_dG']

        # If specified, add the phosphorothioate bonds to the oligo.
        if self.config.num_phosphorothioate_bonds > 0:
            oligo_seq = self.add_phosphorothioate_bonds(
                    oligo_seq, self.config.num_phosphorothioate_bonds)

        # Calculate the predicted replacement efficiency.
        predicted_replacement_efficiency = (
                self.calc_predicted_replacement_efficiency())

        return OligoResult(
            oligo_target.target_id,
            oligo_target.start + 1,
            oligo_target.end + 1,
            oligo_target.strand,
            oligo_target.replichore,
            oligo_target.mutation_type,
            ss_dG,
            self.config.oligo_size,
            oligo_seq,
            oligo_target.original_seq,
            oligo_target.mutation_seq,
            predicted_replacement_efficiency
        )


    def get_candidate_block_seq(self, oligo_target):
        """Returns the sequence block from which we will chose the oligo.

        This takes into account the boundaries set in the config as well
        as the target specs. The method that actually searches for the optimal
        oligo sequence window within this block must stay
        config.oligo_end_buffer_distance away from each end of this
        candidate.

        The polarity of this block sequence is in the direction of the resulting
        oligo.
        """
        # The strategy for determining the candidate block is:
        #     1) Determine upstream and downstream bounds of the block:
        #           Use the [BioPython] convention of positions being relative
        #           to the sense strand.
        #     2) Determine correct polarity:
        #           If the target lies on replichore 1 and (+) strand or on
        #           replichore 2 and on the (-) strand, then the sense strand
        #           is always targeted (and the oligo is on the anti-sense
        #           strand). If the target lies on replichore 1 and (-) or on
        #           replichore 2 and on the (+) strand, then the anti-sense
        #           strand is always targeted (and the oligo lies on the sense
        #           strand).


        ### 1) Determine bounds.

        # The downstream bound (using real numbers to help think about it) is
        # determined by the oligo that would start at the downstream-most
        # position possible, or ~15 base pairs upstream of the mutation region,
        # and end 90 base pairs down stream. We take into account any bases
        # added or deleted at the target:

        downstream_bound = (oligo_target.end + (
                self.config.oligo_size -
                self.config.oligo_end_buffer_distance -
                len(oligo_target.mutation_seq) -
                1))
        downstream_bound = min(downstream_bound,
                len(self.config.genome_record.seq))

        # The upstream bound follows a similar argument.
        upstream_bound = (oligo_target.start - (
                self.config.oligo_size -
                        self.config.oligo_end_buffer_distance -
                        len(oligo_target.mutation_seq)))
        upstream_bound = max(0, upstream_bound)

        # Prepare the mutation sequence to have the correct polarity relative
        # to the block.
        if oligo_target.strand == 1:
            mutation_seq_forward_strand = oligo_target.mutation_seq
        else:
            mutation_seq_forward_strand = (
                    oligo_target.mutation_seq.reverse_complement())

        # Now put together the sense strand version of the block seq. We'll
        # figure out correct polarity next.
        block_seq = (
                self.config.genome_record.seq[
                        upstream_bound:oligo_target.start] +
                mutation_seq_forward_strand +
                self.config.genome_record.seq[oligo_target.end:downstream_bound]
        )


        ### 2) Determine the polarity of the oligo return the appropriate
        ###    block sequence.

        oligo_sense = self.determine_oligo_sense(oligo_target)
        if oligo_sense == 1:
            return block_seq
        else:
            return block_seq.reverse_complement()


    def determine_oligo_sense(self, oligo_target):
        """Determine the sense of the oligo target.

        Args:
            oligo_target: An OligoTarget object.

        Returns:
            Either 1 or -1 to indicate sense that oligo will have.
        """
        return -1 if oligo_target.replichore == 1 else 1


    def determine_oligo_from_block(self, block_seq):
        """Starts centered on the candidate block and wiggles outwards
        until a window is found that satisfies the secondary structure free
        energy constraint.
        """
        initial_midpoint = len(block_seq) / 2
        current_midpoint = initial_midpoint
        wiggle = 1

        # Calculate this once.
        downstream_half_oligo_size = int(
                math.floor(self.config.oligo_size / 2.0))
        upstream_half_oligo_size = int(
                math.ceil(self.config.oligo_size / 2.0))

        # Make the initial cut and set the base for the best oligo.
        upstream_block_cut = (current_midpoint - downstream_half_oligo_size)
        downstream_block_cut = (current_midpoint + upstream_half_oligo_size)
        oligo_seq = block_seq[upstream_block_cut:downstream_block_cut]
        assert self.config.oligo_size == len(oligo_seq), (
                "Expected: %d, Actual: %d" %
                        (self.config.oligo_size, len(oligo_seq)))
        ss_dG = self.get_ss_free_energy(oligo_seq)
        best_oligo_candidate = {
                'oligo_seq': oligo_seq,
                'ss_dG': ss_dG
        }

        # Keep track of the range explored for test/debug.
        debug_midpoint_range_explored = [current_midpoint, current_midpoint]

        while ss_dG < self.config.min_ss_dG:
            current_midpoint = initial_midpoint + wiggle
            upstream_block_cut = (current_midpoint - downstream_half_oligo_size)
            downstream_block_cut = (current_midpoint + upstream_half_oligo_size)

            # Update debug range explored.
            if current_midpoint < debug_midpoint_range_explored[0]:
                debug_midpoint_range_explored[0] = current_midpoint
            elif current_midpoint > debug_midpoint_range_explored[1]:
                debug_midpoint_range_explored[1] = current_midpoint

            # Check if we've hit a boundary of the space to explore.
            # TODO: Right now we're missing potentially checking the other
            # boundary which is okay since we'll get a close-enough answer
            # but can't hurt to fix.
            if upstream_block_cut < 0 or downstream_block_cut > len(block_seq):
                break

            # Make an oligo cut and calculate the free energy.
            oligo_seq = block_seq[upstream_block_cut:downstream_block_cut]
            assert self.config.oligo_size == len(oligo_seq), (
                    "Expected: %d, Actual: %d" %
                            (self.config.oligo_size, len(oligo_seq)))
            ss_dG = self.get_ss_free_energy(oligo_seq)
            if ss_dG > best_oligo_candidate['ss_dG']:
                best_oligo_candidate = {
                    'oligo_seq': oligo_seq,
                    'ss_dG': ss_dG
                }

            # Update the wiggle for the next iteration.
            if wiggle > 0:
                wiggle *= -1
            else:
                wiggle *= -1
                wiggle += 1

        # Return what we have at this point.
        return best_oligo_candidate, debug_midpoint_range_explored


    def get_ss_free_energy(self, seq):
        """Returns the free energy for the given sequence as a float.
        """
        HYBRID_SS_MIN_CMD_ROOT = 'hybrid-ss-min --NA=DNA --energyOnly -q '
        full_cmd = HYBRID_SS_MIN_CMD_ROOT + str(seq)
        p = subprocess.Popen(full_cmd.split(), stdout=subprocess.PIPE)
        return float(p.stdout.read().strip())


    def add_phosphorothioate_bonds(self, oligo_seq, num_thio_bonds):
        """Adds '*' symbols to the resulting oligo, which are recognized by
        IDT as placeholders for phosphorothioate bonds.

        Returns a copy of oligo_seq with the bonds added.
        """
        # Insert a '*' at every other base, starting after the first
        # base in the oligo.
        thio_indeces = [2 * x + 1 for x in range(num_thio_bonds)]

        new_seq_len = len(oligo_seq) + num_thio_bonds

        orig_seq_ptr = 0

        new_seq = ''

        for new_seq_idx in range(new_seq_len):
            if new_seq_idx in thio_indeces:
                new_seq += PHOSPHOROTHIOATE_SYMBOl
            else:
                new_seq += oligo_seq[orig_seq_ptr]
                orig_seq_ptr += 1

        assert new_seq_len == len(new_seq)
        return new_seq


    def calc_predicted_replacement_efficiency(self):
        """Calculates the predicted replacement efficiency.
        """
        # TODO
        return 0.0


class OligoWriter(object):
    """Object that writes a list of oligos to output in the specified format.
    """

    # Header fields for default output format.
    # NOTE: These match the object attributes of OligoResult to make writing more
    # convenient.
    DEFAULT_OUTPUT_OLIGO_FIELD_NAMES = [
            'target_id',
            'start',
            'end',
            'strand',
            'replichore',
            'mutation_type',
            'ss_dG',
            'oligo_size',
            'oligo_seq',
            'original_seq',
            'mutation_seq',
            'predicted_replacement_efficiency'
    ]

    @classmethod
    def write_default(cls, oligo_result_list, oligo_output_file):
        """Write list of result oligos to file.

        Args:
            oligo_result_list: List of OligoResult objects.
            oligo_output_file: Filename or filehandle.
        """
        if isinstance(oligo_output_file, str):
            csvfile = open(oligo_output_file, 'w')
        else:
            csvfile = oligo_output_file

        writer = csv.DictWriter(csvfile, cls.DEFAULT_OUTPUT_OLIGO_FIELD_NAMES)
        writer.writeheader()
        for oligo_result in oligo_result_list:
            writer.writerow(oligo_result.__dict__)


    IDT_OUTPUT_OLIGO_FIELD_NAMES = [
            'target_id',
            'well',
            'sequence'
    ]

    @classmethod
    def write_idt(cls, oligo_result_list, oligo_output_filename):
        """Write list of result oligos to file in IDT format.
        """
        class WellIdGenerator(object):
            """Generates 96-plate well ids from A1 ... A12, ..., H1 ... H12
            """

            LETTER_TRANSITION_TABLE = {
                    'A': 'B',
                    'B': 'C',
                    'C': 'D',
                    'D': 'E',
                    'E': 'F',
                    'F': 'G',
                    'G': 'H',
                    'H': 'A',
            }


            def __init__(self):
                self.letter = 'A'
                self.number = 1


            def __iter__(self):
                return self


            def next(self):
                # Create the current return value.
                current_id = self.letter + "%02d" % (self.number,)

                # Bump the state.
                if self.number == 12:
                    self.letter = self.LETTER_TRANSITION_TABLE[self.letter]

                if self.number == 12:
                    self.number = 1
                else:
                    self.number += 1

                # Return current.
                return current_id


        with open(oligo_output_filename, 'w') as csvfile:
            writer = csv.DictWriter(csvfile,
                    cls.IDT_OUTPUT_OLIGO_FIELD_NAMES)
            writer.writeheader()
            for oligo_result, well_id in zip(oligo_result_list, WellIdGenerator()):
                writer.writerow({
                    'target_id': oligo_result.target_id,
                    'well': well_id,
                    'sequence': oligo_result.oligo_seq
                })


def create_parser():
    """Returns an argparse.ArgumentParser with the options we support.
    """
    # Parse command line arguments.
    # NOTE: These all have default values or are flags.
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref_genome", default=DEFAULT_GENOME,
            help="Manually specified reference genome. "
                    "[Default: %s]" % DEFAULT_GENOME)
    parser.add_argument("--replication_origin", type=tuple,
            default=DEFAULT_REPLICATION_ORIGIN,
            help="Replication origin interval. "
                    "[Default: %s]" % (DEFAULT_REPLICATION_ORIGIN,))
    parser.add_argument("--replication_terminus", type=tuple,
            default=DEFAULT_REPLICATION_TERMINUS,
            help="Replication terminus interval. "
                    "[Default: %s]" % (DEFAULT_REPLICATION_TERMINUS,))
    parser.add_argument("--input_targets",
            default=DEFAULT_OLIGO_TARGET_FILENAME,
            help="Input file containing oligo target specs. "
                    "[Default: %s]" % DEFAULT_OLIGO_TARGET_FILENAME)
    parser.add_argument("--output",
            default=DEFAULT_OLIGO_RESULT_OUTPUT_FILENAME,
            help="Output destination file. "
                    "[Default: %s]" % DEFAULT_OLIGO_RESULT_OUTPUT_FILENAME)
    parser.add_argument("--output_format", default='default',
            choices=['default', 'idt'],
            help="Output format.")
    parser.add_argument("--oligo_size", default=DEFAULT_OLIGO_SIZE,  type=int,
            help="Size of the oligo to create. [Default: %d]" %
                    DEFAULT_OLIGO_SIZE)
    parser.add_argument("--min_ss_dG", default=DEFAULT_MIN_SS_DG, type=float,
            help="Minimum free energy for designed oligos. [Default: %.1f]" %
                    DEFAULT_MIN_SS_DG)
    parser.add_argument("--mut_loc_max", default=DEFAULT_MUT_LOC_MAX, type=int,
            help="Distance between mutation and end of oligo. [Default: %d]" %
                    DEFAULT_MUT_LOC_MAX)
    parser.add_argument("--num_thio", default=DEFAULT_NUM_PHOSPHOROTHIOATE,
            type=int,
            help="Number of phosphorothioate bonds to add. [Default: %d]" %
                    DEFAULT_NUM_PHOSPHOROTHIOATE)
    parser.add_argument("--manually_calc_replichore", action='store_true',
            help="Whether to expect a manually provided replichore for each "
                    "oligo target. Otherwise it is automatically calculated.")
    return parser


def main():
    """Main logic of the script.
    """
    # Parse commandline arguments.
    parser = create_parser()
    args = parser.parse_args()

    # Configure this run from these inputs.
    optMAGE_config = OptMAGEConfig.build_from_args(args)

    # Parse the target input and create list of targets.
    oligo_target_list = OligoTarget.parse_input_and_create_list(
            optMAGE_config, args.input_targets)

    # Generate the oligos.
    oligo_generator = OligoGenerator(optMAGE_config)
    oligo_result_list = [
            oligo_generator.generate_oligo(oligo_target)
            for oligo_target in oligo_target_list]

    # Write the results to the output file.
    if args.output_format == 'idt':
        OligoWriter.write_idt(oligo_result_list, args.output)
    else:
        OligoWriter.write_default(oligo_result_list, args.output)


if __name__ == '__main__':
    main()
