"""
OptMAGE Unit Tests
"""

import os
import unittest

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from optmage.oligo_designer import DEFAULT_GENOME
from optmage.oligo_designer import DEFAULT_OLIGO_SIZE
from optmage.oligo_designer import DEFAULT_MIN_SS_DG
from optmage.oligo_designer import DEFAULT_MUT_LOC_MAX
from optmage.oligo_designer import DEFAULT_NUM_PHOSPHOROTHIOATE
from optmage.oligo_designer import DEFAULT_REPLICATION_ORIGIN
from optmage.oligo_designer import DEFAULT_REPLICATION_TERMINUS
from optmage.oligo_designer import OligoGenerator
from optmage.oligo_designer import OligoTarget
from optmage.oligo_designer import OptMAGEConfig


PWD = os.path.dirname(os.path.realpath(__file__))
TEST_DATA = os.path.join(PWD, 'test_data')
DEFAULT_REF_GENOME = os.path.join(TEST_DATA, 'mg1655.fasta')


class MockObject(object):
    """Generic test object for passing to unit tests.
    """
    pass


class TestOptMAGE(unittest.TestCase):
    """Tests for optMAGE.
    """

    def setUp(self):
        """General setUp routines.
        """
        # Get the reference genome.
        with open(DEFAULT_REF_GENOME) as fh:
            self.ref_genome = SeqIO.read(fh, 'fasta')

        # Create a default config.
        mock_args = MockObject()
        mock_args.oligo_size = DEFAULT_OLIGO_SIZE
        mock_args.min_ss_dG = DEFAULT_MIN_SS_DG
        mock_args.mut_loc_max = DEFAULT_MUT_LOC_MAX
        mock_args.num_thio = DEFAULT_NUM_PHOSPHOROTHIOATE
        mock_args.manually_calc_replichore = False
        mock_args.ref_genome = DEFAULT_GENOME
        mock_args.replication_origin = DEFAULT_REPLICATION_ORIGIN
        mock_args.replication_terminus = DEFAULT_REPLICATION_TERMINUS
        self.config = OptMAGEConfig.build_from_args(mock_args)

    def test_oligo_generator__get_candidate_block_seq(self):
        """Tests the part of the mage generator that determines the candidate
        block seq.
        """
        OLIGO_SIZE = 90
        self.config.oligo_size = OLIGO_SIZE

        OLIGO_END_BUFFER_DISTANCE = 20
        self.config.oligo_end_buffer_distance = OLIGO_END_BUFFER_DISTANCE

        # Set this arbitrarily low so we pretty much guarantee we get the
        # oligo sequence centered around the mutation.
        self.config.min_ss_dG = -200

        # Ignore phosphorothioates for now.
        self.config.num_phosphorothioate_bonds = 0

        oligo_generator = OligoGenerator(self.config)

        # Replichore 1 and negative strand means the anti-sense strand will
        # be targeted, so the oligo will have a positive sense.
        REF = 'C'
        params = {
            'target_id': 'test',
            'replichore': 'NA',
            'strand': '+1',
            'start': 2216229,
            'end': 2216230,
            'mutation_type': 'R'
        }
        oligo_target = OligoTarget(self.config, params)

        # Test getting the candidate block sequence.
        formatted_block_seq = str(
                oligo_generator.get_candidate_block_seq(oligo_target)).upper()
        UPSTREAM_OF_MUT = 'CAACAACCAGCGCCACAGCGGATGCGTGGAGATTCGGCGGATGGCATCGCTACAGGCCAGCAATGCCAG'
        DOWNSTREAM_OF_MUT = 'GCCGCAGCCAGCCAGAAACCACTGCCGAGGCTGGTACGCGCCAGCGCACTGCCATTTTGCGCCAGTTG'
        EXPECTED_BLOCK_SEQ = UPSTREAM_OF_MUT + REF + DOWNSTREAM_OF_MUT
        self.assertEqual(len(EXPECTED_BLOCK_SEQ), len(formatted_block_seq))
        self.assertEqual(EXPECTED_BLOCK_SEQ, formatted_block_seq)

    def test_oligo_generator__determine_oligo_from_block(self):
        """Tests that getting the oligo from block seq works.
        """
        OLIGO_SIZE = 90
        self.config.oligo_size = OLIGO_SIZE

        OLIGO_END_BUFFER_DISTANCE = 20
        self.config.oligo_end_buffer_distance = OLIGO_END_BUFFER_DISTANCE

        # Ensures full space is explored
        self.config.min_ss_dG = 100

        # Ignore phosphorothioates for now.
        self.config.num_phosphorothioate_bonds = 0

        oligo_generator = OligoGenerator(self.config)

        # Use seq cenetered aroud position 2216229.
        REF = 'C'
        UPSTREAM_OF_MUT = 'AACAACCAGCGCCACAGCGGATGCGTGGAGATTCGGCGGATGGCATCGCTACAGGCCAGCAATGCCAG'
        DOWNSTREAM_OF_MUT = 'GCCGCAGCCAGCCAGAAACCACTGCCGAGGCTGGTACGCGCCAGCGCACTGCCATTTTGCGCCAGTTGG'
        BLOCK_SEQ = UPSTREAM_OF_MUT + REF + DOWNSTREAM_OF_MUT

        best_oligo_candidate, midpoint_range_explored = (
                oligo_generator.determine_oligo_from_block(BLOCK_SEQ))

        # Assert that we explored the full space (100 min dG ensures this.)
        # NOTE: We are aware we miss one, but okay with it for simplicity.
        EXPECTED_MIDPOINT_RANGE = [45, len(BLOCK_SEQ) - 44]
        self.assertEqual(EXPECTED_MIDPOINT_RANGE, midpoint_range_explored)


    def test_oligo_generator__from_reference__target_forward_strand(self):
        """Test for synthesizing an oligo that targets the forward strand,
        meaning a reverse-sense oligo.
        """
        OLIGO_SIZE = 90
        self.config.oligo_size = OLIGO_SIZE

        OLIGO_END_BUFFER_DISTANCE = 20
        self.config.oligo_end_buffer_distance = OLIGO_END_BUFFER_DISTANCE

        # Set this arbitrarily low so we pretty much guarantee we get the
        # oligo sequence centered around the mutation.
        self.config.min_ss_dG = -200

        # Ignore phosphorothioates for now.
        self.config.num_phosphorothioate_bonds = 0

        oligo_generator = OligoGenerator(self.config)

        # Replichore 1 and negative strand means the anti-sense strand will
        # be targeted, so the oligo will have a positive sense.
        params = {
            'target_id': 'r_3_set1_der_2635317',
            'replichore': 'NA',
            'strand': -1,
            'start': 2635317,
            'end': 2635318,
            'mutation_type': 'R'
        }
        oligo_target = OligoTarget(self.config, params)

        # Test getting the candidate block sequence.
        formatted_block_seq = str(
                oligo_generator.get_candidate_block_seq(oligo_target)).upper()
        EXPECTED_BLOCK_SEQ = 'CGACCGTACTTACGGTCACGAGTCAGACCCGGGAAATCCGCAACCAGCGCATCTCGGGTGCGAGTTAGACGGTTAAATAACGTGGATTTTCCTACGTTAGGGCGCCCGACAAGCGCGACCACAGGTACCATGTTTAAA'
        self.assertEqual(EXPECTED_BLOCK_SEQ, formatted_block_seq)

        # Test getting the actual oligo seq.
        formatted_oligo_seq = str(
                oligo_generator.generate_oligo(oligo_target).oligo_seq).upper()
        EXPECTED_OLIGO_SEQ = 'AGACCCGGGAAATCCGCAACCAGCGCATCTCGGGTGCGAGTTAGACGGTTAAATAACGTGGATTTTCCTACGTTAGGGCGCCCGACAAGC'
        self.assertEqual(OLIGO_SIZE, len(formatted_oligo_seq))
        self.assertEqual(EXPECTED_OLIGO_SEQ, formatted_oligo_seq)

    def test_oligo_generator__from_reference__target_reverse_strand(self):
        """Test for synthesizing an oligo that targets the reverse strand,
        meaning a forward-sense oligo.
        """
        OLIGO_SIZE = 90
        self.config.oligo_size = OLIGO_SIZE

        OLIGO_END_BUFFER_DISTANCE = 20
        self.config.oligo_end_buffer_distance = OLIGO_END_BUFFER_DISTANCE

        # Set this arbitrarily low so we pretty much guarantee we get the
        # oligo sequence centered around the mutation.
        self.config.min_ss_dG = -200

        # Ignore phosphorothioates for now.
        self.config.num_phosphorothioate_bonds = 0

        oligo_generator = OligoGenerator(self.config)

        # Replichore 1 and positive strand means the anti-sense strand will
        # be targeted, so the oligo will have a positive sense.
        params = {
            'target_id': 'r_1_set1_ftsA_104352',
            'replichore': 'NA',
            'strand': '+1',
            'start': 104352,
            'end': 104353,
            'mutation_type': 'R'
        }
        oligo_target = OligoTarget(self.config, params)

        # Test getting the candidate block sequence.
        formatted_block_seq = str(
                oligo_generator.get_candidate_block_seq(oligo_target)).upper()
        EXPECTED_BLOCK_SEQ = 'CCGGATTCTTGATCCCTTCCTGATAGTCAATCGCATACTCTTGCGGGATCACATGCAGCACACGATGCTCATCGCGCACACGCACCGATTTCGCGGTATGGACGACGTTTTCCACATCTTCTTGCGTCACTTCTTCTT'
        self.assertEqual(EXPECTED_BLOCK_SEQ, formatted_block_seq)

        # Test getting the actual oligo seq.
        formatted_oligo_seq = str(
                oligo_generator.generate_oligo(oligo_target).oligo_seq).upper()
        EXPECTED_OLIGO_SEQ = 'AGTCAATCGCATACTCTTGCGGGATCACATGCAGCACACGATGCTCATCGCGCACACGCACCGATTTCGCGGTATGGACGACGTTTTCCA'
        self.assertEqual(OLIGO_SIZE, len(formatted_oligo_seq))
        self.assertEqual(EXPECTED_OLIGO_SEQ, formatted_oligo_seq)

    def test_determine_oligo_sense(self):
        """Tests for OligoGenerator.determine_oligo_sense().
        """
        OLIGO_GENERATOR = OligoGenerator(self.config)

        # Replichore = 1, strand = +1.
        self.assertEqual(-1, OLIGO_GENERATOR.determine_oligo_sense(
                OligoTarget(self.config, {
                    'target_id': 'r_1_set1_ftsA_104352',
                    'replichore': 'NA',
                    'strand': '+1',
                    'start': 104352,
                    'end': 104353,
                    'mutation_type': 'R'
                })))

        # Replichore = 1, strand = -1.
        self.assertEqual(-1, OLIGO_GENERATOR.determine_oligo_sense(
                OligoTarget(self.config, {
                    'target_id': 'r_8_set1_surA_53597',
                    'replichore': 'NA',
                    'strand': '-1',
                    'start': 53597,
                    'end': 53598,
                    'mutation_type': 'R'
                })))

        # Replichore = 2, strand = +1.
        self.assertEqual(1, OLIGO_GENERATOR.determine_oligo_sense(
                OligoTarget(self.config, {
                    'target_id': 'r_102_set5_csiR_2794168',
                    'replichore': 'NA',
                    'strand': +1,
                    'start': 2794168,
                    'end': 2794169,
                    'mutation_type': 'R'
                })))

        # Replichore = 2, strand = -1.
        self.assertEqual(1, OLIGO_GENERATOR.determine_oligo_sense(
                OligoTarget(self.config, {
                    'target_id': 'r_3_set1_der_2635317',
                    'replichore': 'NA',
                    'strand': -1,
                    'start': 2635317,
                    'end': 2635318,
                    'mutation_type': 'R'
                })))

    def test_free_energy_optmization(self):
        """Tests that the oligo search optimizes the free energy and scans
        both left and right of the mutation midpoint.
        """
        self.config.should_calc_replichore = False
        OLIGO_SIZE = 20
        self.config.oligo_size = OLIGO_SIZE
        OLIGO_END_BUFFER_DISTANCE = 2
        self.config.oligo_end_buffer_distance = OLIGO_END_BUFFER_DISTANCE

        OLIGO_GENERATOR = OligoGenerator(self.config)

        ### Test that the window slides downstream.
        RAW_SEQ_1 = 'GGGGGGGGGGCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
        SEQ_OBJ_1 = Seq(RAW_SEQ_1, generic_dna)
        GENOME_RECORD_1 = SeqRecord(SEQ_OBJ_1)
        self.config.genome_record = GENOME_RECORD_1
        OLIGO_TARGET = OligoTarget(self.config, {
                'target_id': '1',
                'replichore': 2,
                'strand': -1,
                'start': 35,
                'end': 36,
                'mutation_type': 'R'
        })
        oligo_result = OLIGO_GENERATOR.generate_oligo(OLIGO_TARGET)
        self.assertTrue(oligo_result.ss_dG > DEFAULT_MIN_SS_DG)
        EXPECTED_SEQ = 'A*A*AAAAAAAAAAAAAAAAAA'
        self.assertEqual(EXPECTED_SEQ, str(oligo_result.oligo_seq).upper())

        ### Test that the window slides upstream.
        RAW_SEQ_2 = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGCCCCCCCCCCCCC'
        SEQ_OBJ_2 = Seq(RAW_SEQ_2, generic_dna)
        GENOME_RECORD_2 = SeqRecord(SEQ_OBJ_2)
        self.config.genome_record = GENOME_RECORD_2
        OLIGO_TARGET = OligoTarget(self.config, {
                'target_id': '1',
                'replichore': 2,
                'strand': -1,
                'start': 22,
                'end': 23,
                'mutation_type': 'R'
        })
        oligo_result = OLIGO_GENERATOR.generate_oligo(OLIGO_TARGET)
        self.assertTrue(oligo_result.ss_dG > DEFAULT_MIN_SS_DG)
        EXPECTED_SEQ = 'A*A*AAAAAAAAAAAAAAAAAA'
        self.assertEqual(EXPECTED_SEQ, str(oligo_result.oligo_seq).upper())

    def test_mutation__positive_strand(self):
        """Test making a mutation relative to the positive strand.
        """
        self.config.should_calc_replichore = False
        OLIGO_SIZE = 7
        self.config.oligo_size = OLIGO_SIZE
        OLIGO_END_BUFFER_DISTANCE = 0
        self.config.oligo_end_buffer_distance = OLIGO_END_BUFFER_DISTANCE
        self.config.num_phosphorothioate_bonds = 0
        self.config.min_ss_dG = -200 # just want the centered window

        OLIGO_GENERATOR = OligoGenerator(self.config)

        RAW_SEQ_1 = 'CGCTAGCCC'
        SEQ_OBJ_1 = Seq(RAW_SEQ_1, generic_dna)
        GENOME_RECORD_1 = SeqRecord(SEQ_OBJ_1)
        self.config.genome_record = GENOME_RECORD_1

        OLIGO_TARGET = OligoTarget(self.config, {
                'target_id': '1',
                'replichore': 2, # so we get an oligo in the positive sense.
                'strand': 1,
                'start': 4,
                'end': 7,
                'mutation_type': 'M',
                'mutation_seq': 'TAA'
        })
        oligo_result = OLIGO_GENERATOR.generate_oligo(OLIGO_TARGET)
        self.assertEqual(OLIGO_SIZE, len(oligo_result.oligo_seq))
        self.assertEqual('TAG', str(oligo_result.original_seq))
        self.assertEqual('TAA', str(oligo_result.mutation_seq))
        EXPECTED_OLIGO_SEQ = 'GCTAACC'
        self.assertEqual(EXPECTED_OLIGO_SEQ,
                str(oligo_result.oligo_seq).upper())

        # Try similar with oligo size 8.
        OLIGO_SIZE = 8
        self.config.oligo_size = OLIGO_SIZE
        oligo_result = OLIGO_GENERATOR.generate_oligo(OLIGO_TARGET)
        self.assertEqual(OLIGO_SIZE, len(oligo_result.oligo_seq))
        EXPECTED_OLIGO_SEQ = 'CGCTAACC'
        self.assertEqual(EXPECTED_OLIGO_SEQ,
                str(oligo_result.oligo_seq).upper())

    def test_mutation__negative_strand(self):
        """Test making a mutation relative to the negative strand.
        """
        self.config.should_calc_replichore = False
        OLIGO_SIZE = 7
        self.config.oligo_size = OLIGO_SIZE
        OLIGO_END_BUFFER_DISTANCE = 0
        self.config.oligo_end_buffer_distance = OLIGO_END_BUFFER_DISTANCE
        self.config.num_phosphorothioate_bonds = 0
        self.config.min_ss_dG = -200 # just want the centered window

        OLIGO_GENERATOR = OligoGenerator(self.config)

        RAW_SEQ_1 = 'CGCTAGCCC'
        SEQ_OBJ_1 = Seq(RAW_SEQ_1, generic_dna)
        GENOME_RECORD_1 = SeqRecord(SEQ_OBJ_1)
        self.config.genome_record = GENOME_RECORD_1

        OLIGO_TARGET = OligoTarget(self.config, {
                'target_id': '1',
                'replichore': 2, # so we get an oligo in the positive sense.
                'strand': -1,
                'start': 4,
                'end': 7,
                'mutation_type': 'M',
                'mutation_seq': 'TTA'
        })
        oligo_result = OLIGO_GENERATOR.generate_oligo(OLIGO_TARGET)
        self.assertEqual(OLIGO_SIZE, len(oligo_result.oligo_seq))
        EXPECTED_OLIGO_SEQ = 'GCTAACC'
        self.assertEqual(EXPECTED_OLIGO_SEQ,
                str(oligo_result.oligo_seq).upper())

        # Try similar with oligo size 8.
        OLIGO_SIZE = 8
        self.config.oligo_size = OLIGO_SIZE
        oligo_result = OLIGO_GENERATOR.generate_oligo(OLIGO_TARGET)
        self.assertEqual(OLIGO_SIZE, len(oligo_result.oligo_seq))
        EXPECTED_OLIGO_SEQ = 'CGCTAACC'
        self.assertEqual(EXPECTED_OLIGO_SEQ,
                str(oligo_result.oligo_seq).upper())

    def test_deletion(self):
        """Test making a deletion.
        """
        self.config.should_calc_replichore = False
        OLIGO_SIZE = 7
        self.config.oligo_size = OLIGO_SIZE
        OLIGO_END_BUFFER_DISTANCE = 2
        self.config.oligo_end_buffer_distance = OLIGO_END_BUFFER_DISTANCE
        self.config.num_phosphorothioate_bonds = 0
        self.config.min_ss_dG = -200 # just want the centered window

        OLIGO_GENERATOR = OligoGenerator(self.config)

        #                  TCGC AGC
        RAW_SEQ_1 = 'TTTTTTTCGCTAGCCCTTTTTTTTTTTTTTTT'
        SEQ_OBJ_1 = Seq(RAW_SEQ_1, generic_dna)
        GENOME_RECORD_1 = SeqRecord(SEQ_OBJ_1)
        self.config.genome_record = GENOME_RECORD_1

        OLIGO_TARGET = OligoTarget(self.config, {
                'target_id': '1',
                'replichore': 2, # so we get an oligo in the positive sense.
                'strand': 1,
                'start': 11,
                'end': 12,
                'mutation_type': 'D',
        })
        oligo_result = OLIGO_GENERATOR.generate_oligo(OLIGO_TARGET)
        self.assertEqual(OLIGO_SIZE, len(oligo_result.oligo_seq))
        EXPECTED_OLIGO_SEQ = 'TCGCAGC'
        EXPECTED_OLIGO_SEQ_ALTERNATE = 'CGCAGCC'
        self.assertTrue(str(oligo_result.oligo_seq).upper() in
                [EXPECTED_OLIGO_SEQ, EXPECTED_OLIGO_SEQ_ALTERNATE],
                'Got: ' + str(oligo_result.oligo_seq).upper())

        # Try similar with oligo size 8.
        OLIGO_SIZE = 8
        self.config.oligo_size = OLIGO_SIZE
        oligo_result = OLIGO_GENERATOR.generate_oligo(OLIGO_TARGET)
        self.assertEqual(OLIGO_SIZE, len(oligo_result.oligo_seq))
        EXPECTED_OLIGO_SEQ = 'TCGCAGCC'
        EXPECTED_OLIGO_SEQ_ALTERNATE = 'TTCGCAGC'
        self.assertTrue(str(oligo_result.oligo_seq).upper() in
                [EXPECTED_OLIGO_SEQ, EXPECTED_OLIGO_SEQ_ALTERNATE])

        ### Test bigger deletion.
        OLIGO_SIZE = 7
        self.config.oligo_size = OLIGO_SIZE
        OLIGO_TARGET = OligoTarget(self.config, {
                'target_id': '1',
                'replichore': 2, # so we get an oligo in the positive sense.
                'strand': 1,
                'start': 11,
                'end': 14,
                'mutation_type': 'D',
        })
        oligo_result = OLIGO_GENERATOR.generate_oligo(OLIGO_TARGET)
        self.assertEqual(OLIGO_SIZE, len(oligo_result.oligo_seq))
        EXPECTED_OLIGO_SEQ = 'TCGCCCC'
        EXPECTED_OLIGO_SEQ_ALTERNATE = 'TTCGCCC'
        self.assertTrue(str(oligo_result.oligo_seq).upper() in
                [EXPECTED_OLIGO_SEQ, EXPECTED_OLIGO_SEQ_ALTERNATE])

        # Try similar with oligo size 8.
        OLIGO_SIZE = 8
        self.config.oligo_size = OLIGO_SIZE
        oligo_result = OLIGO_GENERATOR.generate_oligo(OLIGO_TARGET)
        self.assertEqual(OLIGO_SIZE, len(oligo_result.oligo_seq))
        EXPECTED_OLIGO_SEQ = 'TTCGCCCC'
        EXPECTED_OLIGO_SEQ_ALTERNATE = 'TCGCCCCT'
        self.assertTrue(str(oligo_result.oligo_seq).upper() in
                [EXPECTED_OLIGO_SEQ, EXPECTED_OLIGO_SEQ_ALTERNATE])

    def test_mutation__full_genome__positive_strand(self):
        self.config.num_phosphorothioate_bonds = 0
        OLIGO_GENERATOR = OligoGenerator(self.config)
        OLIGO_TARGET = OligoTarget(self.config, {
                'target_id': '1',
                'replichore': 2, # so we get an oligo in the positive sense.
                'strand': 1,
                'start': 2216229,
                'end': 2216230,
                'mutation_type': 'M',
                'mutation_seq': 'T'
        })
        oligo_result = OLIGO_GENERATOR.generate_oligo(OLIGO_TARGET)
        self.assertEqual(DEFAULT_OLIGO_SIZE, len(oligo_result.oligo_seq))
        self.assertEqual('C', str(oligo_result.original_seq).upper())
        self.assertEqual('T', str(oligo_result.mutation_seq).upper())

    def test_input_accepts_strings_or_numbers(self):
        """Input might be parsed from file so should handle numbers as
        strings.
        """
        self.config.num_phosphorothioate_bonds = 0
        OLIGO_GENERATOR = OligoGenerator(self.config)
        OLIGO_TARGET = OligoTarget(self.config, {
                'target_id': '1',
                'replichore': 2, # so we get an oligo in the positive sense.
                'strand': 1,
                'start': '2216229', # Testing this.
                'end': 2216230,
                'mutation_type': 'M',
                'mutation_seq': 'T'
        })
        OLIGO_GENERATOR.generate_oligo(OLIGO_TARGET)


if __name__ == '__main__':
    unittest.main()
