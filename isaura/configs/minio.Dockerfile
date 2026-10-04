# MinIO no longer publishes images or binaries, so we build its public AGPL source locally.
FROM golang:1.24-alpine AS build
ARG MINIO_RELEASE=RELEASE.2025-09-07T16-13-09Z
RUN VERSION=$(echo "${MINIO_RELEASE#RELEASE.}" | sed -E 's/T([0-9]+)-([0-9]+)-([0-9]+)Z/T\1:\2:\3Z/') && \
    CGO_ENABLED=0 go install -trimpath -ldflags "-s -w \
      -X github.com/minio/minio/cmd.Version=${VERSION} \
      -X github.com/minio/minio/cmd.ReleaseTag=${MINIO_RELEASE} \
      -X github.com/minio/minio/cmd.CopyrightYear=2025" \
      github.com/minio/minio@${MINIO_RELEASE}

FROM alpine:3.20
RUN apk add --no-cache ca-certificates
COPY --from=build /go/bin/minio /usr/bin/minio
EXPOSE 9000 9001
ENTRYPOINT ["/usr/bin/minio"]
