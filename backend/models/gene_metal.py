import uuid
from datetime import datetime

from sqlalchemy import DateTime, ForeignKey, Text, UniqueConstraint, func
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import Mapped, mapped_column

from backend.database.base import Base


class GeneMetal(Base):
    __tablename__ = "gene_metals"

    __table_args__ = (
        UniqueConstraint(
            "gene_id",
            "metal_id",
            name="uq_gene_metal",
        ),
    )

    id: Mapped[uuid.UUID] = mapped_column(
        UUID(as_uuid=True),
        primary_key=True,
        default=uuid.uuid4,
    )
    gene_id: Mapped[uuid.UUID] = mapped_column(
        UUID(as_uuid=True),
        ForeignKey("genes.id"),
        nullable=False,
    )

    metal_id: Mapped[uuid.UUID] = mapped_column(
        UUID(as_uuid=True),
        ForeignKey("heavy_metals.id"),
        nullable=False,
    )

    evidence_type: Mapped[str | None] = mapped_column(
        Text,
        nullable=True,
    )

    evidence_source: Mapped[str | None] = mapped_column(
        Text,
        nullable=True,
    )

    notes: Mapped[str | None] = mapped_column(
        Text,
        nullable=True,
    )

    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        server_default=func.now(),
        nullable=False,
    )

    def __repr__(self) -> str:
        return f"<GeneMetal(gene_id={self.gene_id}, metal_id={self.metal_id})>"
