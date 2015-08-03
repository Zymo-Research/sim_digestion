import unittest
from sim_digestion import sim_digestion

class TestSimDigestion(unittest.TestCase):
    def setUp(self):
        import tempfile

        repeat_seq = 'ACGCAA'
        total_repeat = 10

        self.r1_file = tempfile.NamedTemporaryFile(delete=False)
        self.r2_file = tempfile.NamedTemporaryFile(delete=False)
        self.size_selection_range = (4, 10)
        self.recog_site = '123'
        self.seq = self.recog_site.join([repeat_seq for _ in range(total_repeat)])
        self.cut_pos = 1
        self.fw_cutting_sites = [len(repeat_seq) + self.cut_pos + (len(repeat_seq) + len(self.recog_site)) * i
            for i in range(total_repeat - 1)]
        self.rev_cutting_sites = [s + len(self.recog_site) - self.cut_pos
            for s in self.fw_cutting_sites]

    def tearDown(self):
        import subprocess

        subprocess.check_call(('rm', self.r1_file.name))
        subprocess.check_call(('rm', self.r2_file.name))

    def test_get_bs_seq(self):
        seq = 'GGGGCTTTCCAGAA'
        result = 'GGGGTTTTTTAGAA'

        self.assertEqual(sim_digestion.get_bs_seq(seq), result)

    def test_get_cutting_sites(self):
        self.assertEqual(
            sim_digestion.get_cutting_sites(self.seq, self.recog_site, self.cut_pos),
            (set(self.fw_cutting_sites), set(self.rev_cutting_sites)),
        )

    def test_get_sequence_cutting_sites(self):
        sites_tuple = sim_digestion.get_cutting_sites_in_parallel(
            self.seq,
            [self.recog_site],
            [self.cut_pos],
            max_nt_per_process=10,
        )
        self.assertEqual(
            sites_tuple,
            (set(self.fw_cutting_sites), set(self.rev_cutting_sites)),
        )

    def test_get_fragments(self):
        fragments = sim_digestion.get_fragments(
            self.seq,
            self.fw_cutting_sites,
            self.size_selection_range,
            is_reverse=False,
        )

    def test_write_fastq_files(self):
        fragments = [
            (10, 20, '+', 'A'*5+'C'*5),
            (35, 135, '-', 'AT'*25+'ACGGA'*10),
        ]
        sim_digestion.write_fastq_files(
            self.r1_file.name,
            self.r2_file.name,
            'chr1',
            fragments,
            read_len=5,
        )

if __name__=='__main__':
    unittest.main()
