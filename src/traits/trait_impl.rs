use crate::{
    bed::{BedDataLine, Chrom, Coordinate},
    bedgraph::BedGraphDataLine,
    traits::ToChromStartEndVal,
};

impl<V: Clone> ToChromStartEndVal<V> for BedDataLine<V> {
    fn to_chrom_start_end_val(
        &self,
    ) -> (Chrom, Coordinate, Coordinate, Option<V>) {
        (self.chrom.clone(), self.start, self.end, self.score.clone())
    }
}

impl<V: Clone> ToChromStartEndVal<V> for BedGraphDataLine<V> {
    fn to_chrom_start_end_val(
        &self,
    ) -> (Chrom, Coordinate, Coordinate, Option<V>) {
        (
            self.chrom.clone(),
            self.start,
            self.end_exclusive,
            Some(self.value.clone()),
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::{bed::BedDataLine, traits::ToChromStartEndVal};

    #[test]
    fn test_bed_data_line_to_chorm_start_end_val() {
        macro_rules! test_data_line {
            ($chrom:expr, $start:expr, $end:expr, $score:expr) => {
                let data_line = BedDataLine {
                    chrom: $chrom.to_string(),
                    start: $start,
                    end: $end,
                    name: None,
                    score: $score,
                    strand: None,
                };

                let (chrom, start, end, val) =
                    data_line.to_chrom_start_end_val();

                assert_eq!(&chrom, $chrom);
                assert_eq!(start, $start);
                assert_eq!(end, $end);
                assert_eq!($score, $score);
            };
        }

        test_data_line!("chr1", 2, 10, None::<i32>);
        test_data_line!("chrX", 200, 201, Some(19i32));
        test_data_line!("chrY", 100, 4000, Some(-2i32));
    }
}
