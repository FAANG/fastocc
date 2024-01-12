FROM python:3.10 as requirements-stage

WORKDIR /tmp

RUN pip install poetry

COPY ./pyproject.toml ./poetry.lock* /tmp/

RUN poetry export -f requirements.txt --output requirements.txt --without-hashes


FROM python:3.10

ENV DEBIAN_FRONTEND=noninteractive
RUN apt update && apt install -y bedtools && rm -rf /var/lib/apt/lists/*

WORKDIR /code

COPY --from=requirements-stage /tmp/requirements.txt /code/requirements.txt

RUN pip install --no-cache-dir --upgrade -r /code/requirements.txt

COPY . /code/

ENTRYPOINT ["python", "-m", "fastocc"]
