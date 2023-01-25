#[derive(Clone, Copy, PartialEq, Eq)]
pub enum CallingModel {
    RawPileup,
    RealignedPileup,
    HaplotypeLikelihood
}

pub fn calling_model_to_str(model: CallingModel) -> &'static str {
    match model {
        CallingModel::RawPileup => "raw-pileup",
        CallingModel::RealignedPileup => "realigned-pileup",
        CallingModel::HaplotypeLikelihood => "haplotype-likelihood",
    }
}

pub fn str_to_calling_model(s: &str) -> Result<CallingModel, &'static str> {
    match s {
        "raw-pileup" => Ok(CallingModel::RawPileup),
        "realigned-pileup" => Ok(CallingModel::RealignedPileup),
        "haplotype-likelihood" => Ok(CallingModel::HaplotypeLikelihood),
        _ => Err("Unknown model string")
    }
}
