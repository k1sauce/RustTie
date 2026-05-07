//! `.bt2` index reader and FM-index backward search.
//!
//! The on-disk format is BowTie 2's; the in-memory representation is ours.
//! See `docs/bt2-format.md` and `format.rs` for details.

pub mod build;
pub mod bwt;
pub mod format;
pub mod reader;
pub mod reference;
pub mod search;

pub use format::EbwtParams;
pub use reader::{Bt2Index, RStart};
pub use reference::{BitPairReference, RefRecord};
pub use search::{RefHit, SaRange, backward_search, exact_hits, resolve_text_pos};
