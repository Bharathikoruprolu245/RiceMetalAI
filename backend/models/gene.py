import uuid
from datetime import datetime

from sqlalchemy import DateTime, ForeignKey, Integer, String, Text, func
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import Mapped, mapped_column

from backend.database.base import Base


class Gene(Base):
    __tablename__ = "genes"

    id: Mapped[uuid.UUID] = mapped_column(
        UUID(as_uuid=True),
        primary_key=True,
        default=uuid.uuid4,
    )

    symbol: Mapped[str] = mapped_column(
        String(100),
        nullable=False,
    )

    full_name: Mapped[str | None] = mapped_column(
        Text,
        nullable=True,
    )

    gene_family_id: Mapped[uuid.UUID] = mapped_column(
        UUID(as_uuid=True),
        ForeignKey("gene_families.id"),
        nullable=False,
    )

    species: Mapped[str] = mapped_column(
        String(150),
        nullable=False,
        default="Oryza sativa",
    )

    ncbi_gene_id: Mapped[str] = mapped_column(
        String(50),
        unique=True,
        nullable=False,
    )

    locus_tag: Mapped[str | None] = mapped_column(
        String(100),
        nullable=True,
    )

    chromosome: Mapped[str | None] = mapped_column(
        String(20),
        nullable=True,
    )

    start_position: Mapped[int | None] = mapped_column(
        Integer,
        nullable=True,
    )

    end_position: Mapped[int | None] = mapped_column(
        Integer,
        nullable=True,
    )

    strand: Mapped[str | None] = mapped_column(
        String(10),
        nullable=True,
    )

    gene_length: Mapped[int | None] = mapped_column(
        Integer,
        nullable=True,
    )

    ncbi_assembly: Mapped[str | None] = mapped_column(
        String(100),
        nullable=True,
    )

    ncbi_chromosome_accession: Mapped[str | None] = mapped_column(
        String(100),
        nullable=True,
    )

    ncbi_url: Mapped[str | None] = mapped_column(
        Text,
        nullable=True,
    )

    description: Mapped[str | None] = mapped_column(
        Text,
        nullable=True,
    )

    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        server_default=func.now(),
        nullable=False,
    )

    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        server_default=func.now(),
        onupdate=func.now(),
        nullable=False,
    )

    def __repr__(self) -> str:
        return f"<Gene(symbol={self.symbol!r}, ncbi_gene_id={self.ncbi_gene_id!r})>"
