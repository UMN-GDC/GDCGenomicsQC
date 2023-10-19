<p align="center">
  <i>  A quality control pipeline for genomics data developed by the Masonic Institute of the Developing Brain at the University of Minnesota.</a>.</i>
  <br/>
</p>

![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
![Shell Script](https://img.shields.io/badge/shell_script-%23121011.svg?style=for-the-badge&logo=gnu-bash&logoColor=white)


This is the source code that runs and all the associated services. A quality control pipeline for genomics data developed by the Masonic Institute of the Developing Brain at the University of Minnesota. The pipeline is built utilizing Plink, Liftover,

# Installation
https://3.basecamp.com/5125554/projects/34247145

Please see the [documentation](https://docs.getoutline.com/s/hosting/) for running your own copy of Outline in a production configuration.

If you have questions or improvements for the docs please create a thread in [GitHub discussions](https://github.com/outline/outline/discussions).

# Development

There is a short guide for [setting up a development environment](https://docs.getoutline.com/s/hosting/doc/local-development-5hEhFRXow7) if you wish to contribute changes, fixes, and improvements to Outline.

## Contributing

Outline is built and maintained by a small team – we'd love your help to fix bugs and add features!

Before submitting a pull request _please_ discuss with the core team by creating or commenting in an issue on [GitHub](https://www.github.com/outline/outline/issues) – we'd also love to hear from you in the [discussions](https://www.github.com/outline/outline/discussions). This way we can ensure that an approach is agreed on before code is written. This will result in a much higher liklihood of your code being accepted.

If you’re looking for ways to get started, here's a list of ways to help us improve Outline:

- [Translation](docs/TRANSLATION.md) into other languages
- Issues with [`good first issue`](https://github.com/outline/outline/labels/good%20first%20issue) label
- Performance improvements, both on server and frontend
- Developer happiness and documentation
- Bugs and other issues listed on GitHub

## Architecture

If you're interested in contributing or learning more about the Outline codebase
please refer to the [architecture document](docs/ARCHITECTURE.md) first for a high level overview of how the application is put together.

## Debugging

In development Outline outputs simple logging to the console, prefixed by categories. In production it outputs JSON logs, these can be easily parsed by your preferred log ingestion pipeline.

HTTP logging is disabled by default, but can be enabled by setting the `DEBUG=http` environment variable.

## Tests

We aim to have sufficient test coverage for critical parts of the application and aren't aiming for 100% unit test coverage. All API endpoints and anything authentication related should be thoroughly tested.

To add new tests, write your tests with [Jest](https://facebook.github.io/jest/) and add a file with `.test.js` extension next to the tested code.

```shell
# To run all tests
make test

# To run backend tests in watch mode
make watch
```

Once the test database is created with `make test` you may individually run
frontend and backend tests directly.

```shell
# To run backend tests
yarn test:server

# To run a specific backend test
yarn test:server myTestFile

# To run frontend tests
yarn test:app
```

## Migrations

Sequelize is used to create and run migrations, for example:

```shell
yarn sequelize migration:generate --name my-migration
yarn sequelize db:migrate
```

Or to run migrations on test database:

```shell
yarn sequelize db:migrate --env test
```

# Activity

![Alt](https://repobeats.axiom.co/api/embed/ff2e4e6918afff1acf9deb72d1ba6b071d586178.svg "Repobeats analytics image")

# License

Outline is [BSL 1.1 licensed](LICENSE).
